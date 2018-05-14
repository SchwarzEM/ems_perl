#!/usr/bin/perl -w

# convg_blastpgp_m-6_to_germane_hits.pl
# 
# Erich Schwarz <emsch@its.caltech.edu>, 10/31/01
#
# Some documentation is given at the *end* of this file.
#  
# Basic goal: given an output from blastpgp in 
#   the format "M=6", extract the subsequences 
#   of homologs that are in the convergent set 
#   and that are specifically matched by the 
#   profile.  This is intended to allow efficient 
#   iterative Blast searches without chimeric
#   proteins muddling the iterations.

$first_infile="";
$e_value="";

# $infile="";
# $outfile="";

unless ($#ARGV == 1) 
{
    print "\n";
    print "  This requires, as input, a convergent blastpgp output\n";
    print "  in flat query-anchored format, no identities and blunt ends;\n";
    print "  it will fail on other file types.\n";
    print "\n";
    print "What file is to be extracted?  ";
    $first_infile = <STDIN>;

    print "\n";
    print "  This also requires, as input, a numerical E value (e.g. 1e-03 or\n";
    print "  1e-04) that will be put into psi-blast or gapped-blast scripts.\n";
    print "\n";
    print "What E value is to be used?  ";
    $e_value = <STDIN>;
}
else
{
    $first_infile = $ARGV[0];
    $e_value = $ARGV[1];
}

chomp($first_infile);
chomp($e_value);

$first_outfile=($first_infile.".extract");

$gi_checklist=($first_infile.".hit_gi_nos");
$rough_gi_checklist=($first_infile.".gi_checklist.1");

# Prepare to generate a script to find which hits will see new hits in one BlastP:

$second_blast_script=($first_infile.".2d_blast.sh");
$rough_second_blast_script=($first_infile.".rough_second_blast_script");

# Prepare to generate a script to do a second wave of 30x psi-blast:

$third_blast_script=($first_infile.".2d_30xpsi-blast.sh");
$rough_third_blast_script=($first_infile.".rough_third_blast_script");

# Prepare to generate a filtering script that will allow the results
#   from *.2d_blast.sh to be used to prune .2d_30xpsi-blast.sh:

$filterscript=($first_infile.".filter.sh");  
$rough_filterscript=($first_infile.".rough_filterscript");

# Prepare to generate a script to scan the Ensembl predicted human peptides:

$fourth_blast_script=($first_infile.".ensembl_blast.sh");
$rough_fourth_blast_script=($first_infile.".rough_fourth_blast_script");

# Prepare to generate a script to scan the Ensembl predicted mouse peptides:

$fifth_blast_script=($first_infile.".ens_mouse_blast.sh");
$rough_fifth_blast_script=($first_infile.".rough_fifth_blast_script");


print "\n";
print "The file to be extracted is:         $first_infile\n";
print "Sorted gi numbers of hits will be:   $gi_checklist\n";
print "The extract-file to be generated is: $first_outfile\n";
print "\n";

open (FIRST_INFILE, "$first_infile") || die "Can't open $first_infile. $! \n";
open (FIRST_OUTFILE, ">$first_outfile") || die "Can't open $first_outfile. $! \n";
open (ROUGH_GI_CHECKLIST, ">$rough_gi_checklist") || die "Can't open $gi_checklist. $! \n";

$reading_gi_nos=0; # i.e., not reading gi numbers yet.
$converged_yet=0;
$ignoring_matches=0;
$i=0;
@significant_gi_numbers=();

while (<FIRST_INFILE>) 
{
    if ($_ =~ /Sequences producing significant alignments/) 
    {
        @significant_gi_numbers=();
        $reading_gi_nos=1;
        $reading_matches=0;
        $ignoring_matches=0;
        print "Significant gi numbers, tentatively record.\n";
    }
    elsif ($_ =~ /^gi\|(\d+)\|/ && $reading_gi_nos==1)
    {
        push (@significant_gi_numbers,$1); # push $1 into @significant_gi_numbers;
    }
    elsif ($_ =~ /Sequences not found previously or not previously below threshold/)
    {
        $reading_gi_nos=0;
        print "INsignificant gi numbers, definitely ignore.\n";
    }
    elsif ($_ =~ /^QUERY/ && $converged_yet==0 && $ignoring_matches==0)
    {
            # if no "CONVERGED!" then forget all previous gi nos.
            $reading_gi_nos=0;
            $ignoring_matches=1;
            print "Not converged yet, so ignore all previous gi numbers.\n";
            print "\n";
    }
    elsif ($_ =~ /CONVERGED!/) 
    {
        $converged_yet=1;
        $reading_gi_nos=0;
        print "OK, now I can actually record those significant gi numbers.\n";
        print "\n";
        until ($i == ($#significant_gi_numbers+1)) 
        {
            print ROUGH_GI_CHECKLIST $significant_gi_numbers[$i];
            print ROUGH_GI_CHECKLIST "\n";
            $i = $i+1;
        } 
    }
    elsif ($converged_yet==1 && $_ =~ /^QUERY/) 
    {
        $reading_matches=1;
    }
    elsif ($converged_yet==1 && $reading_matches==1) 
    {
        if ($_ =~ /^\d+\s+\d+\s+.*\w+.*\d+/) 
        { 
            $match_line=$_;
            print FIRST_OUTFILE $match_line;
        }
    }    
}

close FIRST_INFILE;
close FIRST_OUTFILE;
close ROUGH_GI_CHECKLIST;
system `sort $rough_gi_checklist > $gi_checklist`;
system `rm $rough_gi_checklist`;


# Generate a script to find which hits will see new hits in one BlastP:
# Script-specific lines:

open (BLAST_SCRIPT, ">$second_blast_script") 
    || die "Can't open $second_blast_script. $! \n";
open (ROUGH_BLAST_SCRIPT, ">$rough_second_blast_script") 
    || die "Can't open $rough_second_blast_script. $! \n";

open (FILTER_SCRIPT, ">$filterscript") 
    || die "Can't open $filterscript. $! \n";
open (ROUGH_FILTER_SCRIPT, ">$rough_filterscript") 
    || die "Can't open $rough_filterscript. $! \n";

# Somewhat generic code follows:

print BLAST_SCRIPT "#!/bin/bash\n";
print BLAST_SCRIPT "# \n";
print BLAST_SCRIPT "# This script was generated by extracting relevant blastpgp hits \n";
print BLAST_SCRIPT "#     from a converged blastpgp output -- $first_infile --\n";
print BLAST_SCRIPT "#     in \"-m 6\" format.  It is designed to allow each hit to be\n";
print BLAST_SCRIPT "#     tested once for the ability to detect divergent new hits.\n";
print BLAST_SCRIPT "# \n";

print FILTER_SCRIPT "#!/bin/bash\n";
print FILTER_SCRIPT "# \n";
print FILTER_SCRIPT "# This script was generated by extracting relevant blastpgp hits \n";
print FILTER_SCRIPT "#     from a converged blastpgp output -- $first_infile --\n";
print FILTER_SCRIPT "#     in \"-m 6\" format.  It is designed to allow each hit's gapped\n";
print FILTER_SCRIPT "#     BlastP output to be tested for divergent new hits.\n";
print FILTER_SCRIPT "# \n";

$new_residues_of_hit="";
$residues_of_hit="";
@residue_numbers_of_hit=();
$i=0;
until ($i == ($#significant_gi_numbers+1))  
{
    $new_residues_of_hit="";  
    $residues_of_hit="";
    @residue_numbers_of_hit=();
    $fasta_gi = $significant_gi_numbers[$i];
    open (SECOND_INFILE, "$first_outfile") || die "Can't open $first_outfile. $! \n";
    while (<SECOND_INFILE>) 
    {
        if ($_ =~ /^$fasta_gi\s+(\d+)\s+(.*\w+.*)\s+(\d+)/)
        {
            push (@residue_numbers_of_hit,$1);
            push (@residue_numbers_of_hit,$3);
            $new_residues_of_hit=$2;
            $new_residues_of_hit =~ s/\W//g;
            $residues_of_hit = $residues_of_hit . $new_residues_of_hit;
        }
    }
    $i = $i+1;
    $fasta_gi_segment="$fasta_gi"."_aa_"."$residue_numbers_of_hit[0]"."-"."$residue_numbers_of_hit[$#residue_numbers_of_hit]";

    open (SECOND_OUTFILE, ">$fasta_gi_segment") || die "Can't open $fasta_gi_segment. $! \n";

    print SECOND_OUTFILE ">$fasta_gi_segment\n";

    print ROUGH_BLAST_SCRIPT "nice ../blastpgp ";	# Note: this script is designed to be run
								#    from a subdirectory, to avoid clutter.
    print ROUGH_FILTER_SCRIPT "../test_for_fruitful_hits.pl ";

    print ROUGH_BLAST_SCRIPT "-d ../nr ";	# Database [String] -- default = nr
    print ROUGH_BLAST_SCRIPT "-F T ";		# -F  Filter query sequence with SEG? (T/F)
    print ROUGH_BLAST_SCRIPT "-a 2 ";		# -a  Number of processors to use (1 or more)
    print ROUGH_BLAST_SCRIPT "-I T ";		# -I  Show GI's in deflines? (T/F)
    print ROUGH_BLAST_SCRIPT "-e $e_value ";	# -e  Expectation value (E) [default 10.0]
    print ROUGH_BLAST_SCRIPT "-j 1 ";		# -j  Maximum number of passes to use in multipass version (1 or more)
    print ROUGH_BLAST_SCRIPT "-v 1000 ";		# -v  Number of database sequences to show one-line descriptions for
    print ROUGH_BLAST_SCRIPT "-b 1000 ";		# -b  Number of database sequence to show alignments for
    print ROUGH_BLAST_SCRIPT "-m 6 ";		# -m  aln. view 6 = flat query-anchored, no identities and blunt ends
    print ROUGH_BLAST_SCRIPT "-i $fasta_gi_segment ";	# -i  Query File
    print ROUGH_BLAST_SCRIPT "-o $fasta_gi_segment";	# -o  Output File for Alignment
    print ROUGH_BLAST_SCRIPT ".psi-blast-x1.nr_";
    print ROUGH_BLAST_SCRIPT "$e_value";
    print ROUGH_BLAST_SCRIPT ".align-6 ";
    print ROUGH_BLAST_SCRIPT ";\n";


    print ROUGH_FILTER_SCRIPT "$fasta_gi_segment";
    print ROUGH_FILTER_SCRIPT ".psi-blast-x1.nr_";
    print ROUGH_FILTER_SCRIPT "$e_value";
    print ROUGH_FILTER_SCRIPT ".align-6 ";
    print ROUGH_FILTER_SCRIPT "$gi_checklist ";
    print ROUGH_FILTER_SCRIPT "$fasta_gi ";
    print ROUGH_FILTER_SCRIPT "$fasta_gi_segment\n";

#   Note: this defines test_for_fruitful_hits.pl as a script that takes the following arguments:
#       argument 0  =>  m=6 psi-blast 1x output to scan.
#       argument 1  =>  list of pre-detected gi numbers to check against.
#       argument 2  =>  gi number to mark as "having fruitful hits" or not.
#       argument 3  =>  peptide FASTA file to mark as "having fruitful hits" or not.


    while ($residues_of_hit =~ /^(\w{1,60})(.*)/) 
    {
        $residues_of_hit_printline=$1;
        $residues_of_hit_waiting=$2;
        print SECOND_OUTFILE "$residues_of_hit_printline\n";
        $residues_of_hit=$residues_of_hit_waiting;
    }
}

close SECOND_INFILE;
close SECOND_OUTFILE;
close BLAST_SCRIPT;
close ROUGH_BLAST_SCRIPT;
close FILTER_SCRIPT;
close ROUGH_FILTER_SCRIPT;

# end somewhat generic code

# a few more script-specific lines:

system `sort $rough_second_blast_script >> $second_blast_script`;
system `rm $rough_second_blast_script`;
system `chmod +x $second_blast_script`;

system `sort $rough_filterscript >> $filterscript`;
# system `rm $rough_filterscript`;
system `chmod +x $filterscript`;


# Generate a script to do a second wave of 30x psi-blast:
# Script-specific lines:

open (BLAST_SCRIPT, ">$third_blast_script") 
    || die "Can't open $third_blast_script. $! \n";
open (ROUGH_BLAST_SCRIPT, ">$rough_third_blast_script") 
    || die "Can't open $rough_third_blast_script. $! \n";

# Somewhat generic code follows:

print BLAST_SCRIPT "#!/bin/bash\n";
print BLAST_SCRIPT "# \n";
print BLAST_SCRIPT "# This script was generated by extracting relevant blastpgp hits \n";
print BLAST_SCRIPT "#     from a converged blastpgp output -- $first_infile --\n";
print BLAST_SCRIPT "#     in \"-m 6\" format.  It is designed to allow each hit to be\n";
print BLAST_SCRIPT "#     itself run through 30 rounds of psi-blast at a stringency\n";
print BLAST_SCRIPT "#     of $e_value.  IT SHOULD BE HEAVILY FILTERED, to select only for\n";
print BLAST_SCRIPT "#     searches that will in fact give new results on the first round.\n";
print BLAST_SCRIPT "# \n";

$new_residues_of_hit="";
$residues_of_hit="";
@residue_numbers_of_hit=();
$i=0;
until ($i == ($#significant_gi_numbers+1))  
{
    $new_residues_of_hit="";  
    $residues_of_hit="";
    @residue_numbers_of_hit=();
    $fasta_gi = $significant_gi_numbers[$i];
    open (SECOND_INFILE, "$first_outfile") || die "Can't open $first_outfile. $! \n";
    while (<SECOND_INFILE>) 
    {
        if ($_ =~ /^$fasta_gi\s+(\d+)\s+(.*\w+.*)\s+(\d+)/)
        {
            push (@residue_numbers_of_hit,$1);
            push (@residue_numbers_of_hit,$3);
            $new_residues_of_hit=$2;
            $new_residues_of_hit =~ s/\W//g;
            $residues_of_hit = $residues_of_hit . $new_residues_of_hit;
        }
    }
    $i = $i+1;
    $fasta_gi_segment="$fasta_gi"."_aa_"."$residue_numbers_of_hit[0]"."-"."$residue_numbers_of_hit[$#residue_numbers_of_hit]";

    print ROUGH_BLAST_SCRIPT "nice ../blastpgp ";	# Note: this script is designed to be run
								#    from a subdirectory, to avoid clutter.
    print ROUGH_BLAST_SCRIPT "-d ../nr ";	# Database [String] -- default = nr
    print ROUGH_BLAST_SCRIPT "-F T ";		# -F  Filter query sequence with SEG? (T/F)
    print ROUGH_BLAST_SCRIPT "-a 2 ";		# -a  Number of processors to use (1 or more)
    print ROUGH_BLAST_SCRIPT "-I T ";		# -I  Show GI's in deflines? (T/F)
    print ROUGH_BLAST_SCRIPT "-e $e_value ";	# -e  Expectation value (E) [default 10.0]
    print ROUGH_BLAST_SCRIPT "-j 30 ";		# -j  Maximum number of passes to use in multipass version (1 or more)
    print ROUGH_BLAST_SCRIPT "-v 1000 ";	# -v  Number of database sequences to show one-line descriptions for
    print ROUGH_BLAST_SCRIPT "-b 1000 ";	# -b  Number of database sequence to show alignments for
    print ROUGH_BLAST_SCRIPT "-m 6 ";		# -m  aln. view 6 = flat query-anchored, no identities and blunt ends
    print ROUGH_BLAST_SCRIPT "-i $fasta_gi_segment ";	# -i  Query File

    print ROUGH_BLAST_SCRIPT "-o $fasta_gi_segment";	# -o  Output File for Alignment
    print ROUGH_BLAST_SCRIPT ".psi-blast-x30.nr_";
    print ROUGH_BLAST_SCRIPT "$e_value";
    print ROUGH_BLAST_SCRIPT ".align-6 ";

    print ROUGH_BLAST_SCRIPT ";\n";
}

close SECOND_INFILE;
close BLAST_SCRIPT;
close ROUGH_BLAST_SCRIPT;

# end somewhat generic code

# a few more script-specific lines:

system `sort $rough_third_blast_script >> $third_blast_script`;
system `rm $rough_third_blast_script`;
system `chmod +x $third_blast_script`;


# Script-specific lines, to generate a script to scan the Ensembl predicted human peptides: 

open (BLAST_SCRIPT, ">$fourth_blast_script") 
    || die "Can't open $fourth_blast_script. $! \n";
open (ROUGH_BLAST_SCRIPT, ">$rough_fourth_blast_script") 
    || die "Can't open $rough_fourth_blast_script. $! \n";

# more somewhat generic code:

print BLAST_SCRIPT "#!/bin/bash\n";
print BLAST_SCRIPT "# \n";
print BLAST_SCRIPT "# This script was generated by extracting relevant blastpgp hits \n";
print BLAST_SCRIPT "#     from a converged blastpgp output -- $first_infile --\n";
print BLAST_SCRIPT "#     in \"-m 6\" format.  It is designed to allow each hit to be\n";
print BLAST_SCRIPT "#     used to scan the Ensembl predicted human peptide database.\n";
print BLAST_SCRIPT "# \n";

$new_residues_of_hit="";
$residues_of_hit="";
@residue_numbers_of_hit=();
$i=0;
until ($i == ($#significant_gi_numbers+1))  
{
    $new_residues_of_hit="";  
    $residues_of_hit="";
    @residue_numbers_of_hit=();
    $fasta_gi = $significant_gi_numbers[$i];
    open (SECOND_INFILE, "$first_outfile") || die "Can't open $first_outfile. $! \n";
    while (<SECOND_INFILE>) 
    {
        if ($_ =~ /^$fasta_gi\s+(\d+)\s+(.*\w+.*)\s+(\d+)/)
        {
            push (@residue_numbers_of_hit,$1);
            push (@residue_numbers_of_hit,$3);
            $new_residues_of_hit=$2;
            $new_residues_of_hit =~ s/\W//g;
            $residues_of_hit = $residues_of_hit . $new_residues_of_hit;
        }
    }
    $i = $i+1;
    $fasta_gi_segment="$fasta_gi"."_aa_"."$residue_numbers_of_hit[0]"."-"."$residue_numbers_of_hit[$#residue_numbers_of_hit]";

    print ROUGH_BLAST_SCRIPT "nice ../blastpgp ";
    print ROUGH_BLAST_SCRIPT "-d ../ensembl.genscan.fa ";		# Database [String] -- default = nr
    print ROUGH_BLAST_SCRIPT "-F T ";		# -F  Filter query sequence with SEG? (T/F)
    print ROUGH_BLAST_SCRIPT "-a 2 ";		# -a  Number of processors to use (1 or more)
    print ROUGH_BLAST_SCRIPT "-I T ";		# -I  Show GI's in deflines? (T/F)
    print ROUGH_BLAST_SCRIPT "-e 5e-02 ";	# -e  Expectation value (E) [default 10.0]
    print ROUGH_BLAST_SCRIPT "-j 1 ";		# -j  Maximum number of passes to use in multipass version (1 or more)
    print ROUGH_BLAST_SCRIPT "-v 1000 ";		# -v  Number of database sequences to show one-line descriptions for
    print ROUGH_BLAST_SCRIPT "-b 1000 ";		# -b  Number of database sequence to show alignments for
    print ROUGH_BLAST_SCRIPT "-m 6 ";		# -m  aln. view 6 = flat query-anchored, no identities and blunt ends
    print ROUGH_BLAST_SCRIPT "-i $fasta_gi_segment ";	# -i  Query File
    print ROUGH_BLAST_SCRIPT "-o $fasta_gi_segment";	# -o  Output File for Alignment
        print ROUGH_BLAST_SCRIPT ".psi-blast-x1.ensembl.genscan.fa_5e-02.align-6 ";
    print ROUGH_BLAST_SCRIPT ";\n";
}

close SECOND_INFILE;
close BLAST_SCRIPT;
close ROUGH_BLAST_SCRIPT;

# end somewhat-generic code

# more script-specific lines:

system `sort $rough_fourth_blast_script >> $fourth_blast_script`;
system `rm $rough_fourth_blast_script`;
system `chmod +x $fourth_blast_script`;


# Generate a script to scan the Ensembl predicted mouse peptides:
# Script-specific lines:

open (BLAST_SCRIPT, ">$fifth_blast_script") 
    || die "Can't open $fifth_blast_script. $! \n";
open (ROUGH_BLAST_SCRIPT, ">$rough_fifth_blast_script") 
    || die "Can't open $rough_fifth_blast_script. $! \n";

# Somewhat generic code follows:

print BLAST_SCRIPT "#!/bin/bash\n";
print BLAST_SCRIPT "# \n";
print BLAST_SCRIPT "# This script was generated by extracting relevant blastpgp hits \n";
print BLAST_SCRIPT "#     from a converged blastpgp output -- $first_infile --\n";
print BLAST_SCRIPT "#     in \"-m 6\" format.  It is designed to allow each hit to be\n";
print BLAST_SCRIPT "#     used to scan the Ensembl predicted MOUSE peptide database.\n";
print BLAST_SCRIPT "# \n";

$new_residues_of_hit="";
$residues_of_hit="";
@residue_numbers_of_hit=();
$i=0;
until ($i == ($#significant_gi_numbers+1))  
{
    $new_residues_of_hit="";  
    $residues_of_hit="";
    @residue_numbers_of_hit=();
    $fasta_gi = $significant_gi_numbers[$i];
    open (SECOND_INFILE, "$first_outfile") || die "Can't open $first_outfile. $! \n";
    while (<SECOND_INFILE>) 
    {
        if ($_ =~ /^$fasta_gi\s+(\d+)\s+(.*\w+.*)\s+(\d+)/)
        {
            push (@residue_numbers_of_hit,$1);
            push (@residue_numbers_of_hit,$3);
            $new_residues_of_hit=$2;
            $new_residues_of_hit =~ s/\W//g;
            $residues_of_hit = $residues_of_hit . $new_residues_of_hit;
        }
    }
    $i = $i+1;
    $fasta_gi_segment="$fasta_gi"."_aa_"."$residue_numbers_of_hit[0]"."-"."$residue_numbers_of_hit[$#residue_numbers_of_hit]";

    print ROUGH_BLAST_SCRIPT "nice ../blastpgp ";	# Note: this script is designed to be run
								#    from a subdirectory, to avoid clutter.
    print ROUGH_BLAST_SCRIPT "-d ../mouse_genscan.pep ";	# Database [String] -- default = nr
    print ROUGH_BLAST_SCRIPT "-F T ";		# -F  Filter query sequence with SEG? (T/F)
    print ROUGH_BLAST_SCRIPT "-a 2 ";		# -a  Number of processors to use (1 or more)
    print ROUGH_BLAST_SCRIPT "-I T ";		# -I  Show GI's in deflines? (T/F)
    print ROUGH_BLAST_SCRIPT "-e 5e-02 ";	# -e  Expectation value (E) [default 10.0]
    print ROUGH_BLAST_SCRIPT "-j 1 ";		# -j  Maximum number of passes to use in multipass version (1 or more)
    print ROUGH_BLAST_SCRIPT "-v 1000 ";		# -v  Number of database sequences to show one-line descriptions for
    print ROUGH_BLAST_SCRIPT "-b 1000 ";		# -b  Number of database sequence to show alignments for
    print ROUGH_BLAST_SCRIPT "-m 6 ";		# -m  aln. view 6 = flat query-anchored, no identities and blunt ends
    print ROUGH_BLAST_SCRIPT "-i $fasta_gi_segment ";	# -i  Query File
    print ROUGH_BLAST_SCRIPT "-o $fasta_gi_segment";	# -o  Output File for Alignment
        print ROUGH_BLAST_SCRIPT ".psi-blast-x1.mouse_genscan.pep_5e-02.align-6 ";
    print ROUGH_BLAST_SCRIPT ";\n";
}

close SECOND_INFILE;
close BLAST_SCRIPT;
close ROUGH_BLAST_SCRIPT;

# end somewhat generic code

# a few more script-specific lines:

system `sort $rough_fifth_blast_script >> $fifth_blast_script`;
system `rm $rough_fifth_blast_script`;
system `chmod +x $fifth_blast_script`;

# Off into the sunset.

print "\n";
print "My work is done!\n";
print "\n";
