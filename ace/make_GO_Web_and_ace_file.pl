#!/usr/bin/perl

# make_GO_Web_and_ace_file.pl
# 
# Erich Schwarz, 1/28/02
# 
# Purpose: make .ace and GO-format files of GO terms from 
# Tablemaker outputs.  For details, see "WORDY COMMENTS" below.

# #################### run with standard defaults ####################

$defaults = "yes";

print "\n";
print "Run with standard defaults?  [Default is \"yes\"]: ";
$defaults = <STDIN>;
$defaults =~ tr/A-Z/a-z/;
unless ($defaults =~ /no/) 
{
    print "\n";
    print "OK, running with standard defaults...\n";
    $author_address = 'Erich Schwarz <emsch@its.caltech.edu>';
    $seq_list = "sequence.list";
    $pheno_go_list = "pheno2go.list";
    $rnai2go_file = "rnai2go.tab";
    $worm2go_infile = "wormpep2go";
    $gene_assoc_outfile = "gene_association.wb";
    $go2seq_ace = "GO_terms_to_Sequences.ace";
    $onto_hash_file = "ontology.hash";
    $cds_to_locus_file = "wormpep.CDS-to-locus.hash";
    $problem_name_file = "problematic_sequence_names";
    $shoddy_input_file = "shoddy_input_lines";
    $rnai_seqlist = "rnai_specific_sequence_list";
}

# #################### run with hand-entered defaults ####################

elsif ($defaults =~ /no/) 
{
    print "\n";
    print "Please give your name and e-mail address.\n";
    print "Default is \"Erich Schwarz <emsch\@its.caltech.edu>\": ";
    $author_address = <STDIN>;
    chomp($author_address);
    unless ($author_address =~ /\w+/)
    {
        $author_address = 'Erich Schwarz <emsch@its.caltech.edu>';
    }

    print "\n";
    print "Please name the phenotype-to-GO_term list to use (one Phenotype per line).\n";
    print "The default is \"pheno2go.list\": ";
    $pheno_go_list=<STDIN>;
    chomp($pheno_go_list);
    unless ($pheno_go_list =~ /\w+/)
    {
        $pheno_go_list = "pheno2go.list";
    }

    print "\n";
    print "Please name the RNAi -> GO table file to use (one RNAi_result,\n";
    print "   Phenotype, Sequence, and [Protein] per line).\n";
    print "   The default is \"rnai2go.tab\": ";
    $rnai2go_file = <STDIN>;
    chomp($rnai2go_file);
    unless ($rnai2go_file =~ /\w+/)
    {
        $rnai2go_file = "rnai2go.tab";
    }

    print "\n";
    print "Please name the wormpep2go file to use (the default is \"wormpep2go\"): ";
    $wormpep2go_file = <STDIN>;
    chomp ($wormpep2go_file);
    unless ($wormpep2go_file =~ /\w+/)
    {
        $worm2go_infile = "wormpep2go";
    }

    print "\n";
    print "Please name the ontology.hash file to use (default is \"ontology.hash\"): ";
    $comp_file = <STDIN>;
    chomp ($comp_file);
    unless ($comp_file =~ /\w+/)
    {
        $onto_hash_file = "ontology.hash";
    }

    print "\n";
    print "Please name the wormpep.CDS-to-locus.hash to use (default is \"wormpep.CDS-to-locus.hash\"): ";
    $cds_to_locus_file = <STDIN>;
    chomp ($cds_to_locus_file);
    unless ($cds_to_locus_file =~ /\w+/) 
    {
        $cds_to_locus_file = "wormpep.CDS-to-locus.hash";
    }

    print "\n";
    print "Please give the name that the sequence list should have (one line per sequence).\n";
    print "The default is \"sequence.list\": ";
    $seq_list=<STDIN>;
    chomp($seq_list);
    unless ($seq_list =~ /\w+/)
    {
        $seq_list = "sequence.list";
    }

    print "\n";
    print "Please give the name that the gene association output file should have;\n";
    print "The default choice would be \"gene_association.wb\": ";
    $gene_assoc_outfile = <STDIN>;
    chomp($gene_assoc_outfile);
    unless ($gene_assoc_outfile =~ /\w+/)
    {
        $gene_assoc_outfile = "gene_association.wb";
    }

    print "\n";
    print "Please give the name that the GO-term-for-Sequence .ace file should have;\n";
    print "The default choice would be \"GO_terms_to_Sequences.ace\": ";
    $go2seq_ace = <STDIN>;
    chomp($go2seq_ace);
    unless ($go2seq_ace =~ /\w+/)
    {
        $go2seq_ace = "GO_terms_to_Sequences.ace";
    }

    print "\n";
    print "Please give the name that the file for problematic sequence names\n";
    print "    should have (the default would be \"problematic_sequence_names\"): ";
    $problem_name_file = <STDIN>;
    chomp ($problem_name_file);  
    unless ($problem_name_file =~ /\w+/)
    {
        $problem_name_file = "problematic_sequence_names";
    }

    print "\n";
    print "Please give the name that the file for problematic input lines\n";
    print "    should have (the default is \"shoddy_input_lines\"): ";
    $shoddy_input_file = <STDIN>;
    chomp ($shoddy_input_file);
    unless ($shoddy_input_file =~ /\w+/)
    {
        $shoddy_input_file = "shoddy_input_lines";
    }

    print "\n";
    print "Please give the name that the file for RNAi-specific sequences\n";
    print "    should have (the default is \"rnai_specific_sequence_list\"): ";
    $rnai_seqlist = <STDIN>;
    chomp ($rnai_seqlist);
    unless ($rnai_seqlist =~ /\w+/)
    {
        $rnai_seqlist = "rnai_specific_sequence_list";
    }
}

# #################### abort run if defaults not entered ####################

else 
{
    print "Hm.  I didn't get \"yes\" or \"no\" for the defaults.  Try again?\n";
    die;
}

# #################### announce working defaults ####################

print "\n";
print "The author e-mail address is:                $author_address\n";
print "The phenotype-to-GO list is:                 $pheno_go_list\n";
print "The RNAi-to-GO list is:                      $rnai2go_file\n";
print "The wormpep-to-GO file is:                   $worm2go_infile\n";
print "The gene association output file is:         $gene_assoc_outfile\n";
print "The ontology hash file is:                   $onto_hash_file\n";
print "The wormpep.CDS-to-locus.hash is:            $cds_to_locus_file\n";
print "\n";
print "The problematic sequence-name file will be:  $problem_name_file\n";
print "The shoddy input line file will be:          $shoddy_input_file\n";
print "The RNAi-specific seq. list will be:         $rnai_seqlist\n";
print "The sequence list will be:                   $seq_list\n";
print "\n";
print "The gene association output file will be:    $gene_assoc_outfile\n";
print "\n";

# #################### make the date-stamp for GO files ####################

$datestamp = `date -I`;
chomp ($datestamp);
$orig_datestamp = $datestamp;
$datestamp =~ s/-//g;

# #################### generate list of sequences to curate ####################

open (WORM2GO, "$worm2go_infile") || die "Wormpep2go file $worm2go_infile not found. $!\n";
open (ROUGHSEQLIST, ">$seq_list.1") || die "Rough sequence list file $seq_list.1 could not be opened. $!\n";

# format of WORM2GO is:
# AH6.12	WP:CE01453	INTERPRO:IPR000344	GO:0003677	

while (<WORM2GO>) 
{
    if ($_ =~ /^(\w+\.\w+)\W+/) 
    {
    $orfname = $1;
    $orfname =~ tr/a-z/A-Z/;
    print ROUGHSEQLIST "$orfname\n";
    }
}
close WORM2GO;

open (RNA2GO, "$rnai2go_file") || die "RNAi-to-GO file $rnai2go_file not opened. $!\n";
open (ROUGH_RNASEQLIST, ">$rnai_seqlist.1") 
    || die "Rough RNAi-specific sequence list $rnai_seqlist.1 not opened. $!\n";

# format of RNA2GO is:
# "JA:B0019.1"	"WT"	"B0019.1"	"WP:CE19654"

while (<RNA2GO>) 
{
    $inputline = $_;
    $inputline =~ s/\"//g;
    if ($inputline =~ /^\w+\:\w+.+\W+(\w+\.\w+)\W+/) 
    {
    $orfname = $1;
    $orfname =~ tr/a-z/A-Z/;
    print ROUGHSEQLIST "$orfname\n";
    print ROUGH_RNASEQLIST "$orfname\n";
    }
}
close RNA2GO;
close ROUGHSEQLIST;
close ROUGH_RNASEQLIST;

system "sort $seq_list.1 | uniq - > $seq_list" || die "Sequence list $seq_list not created. $!\n";
system "rm $seq_list.1";
system "sort $rnai_seqlist.1 | uniq - > $rnai_seqlist" 
    || die "RNA-specific sequence list $rnai_seqlist not created. $!\n";
print "Sequence lists $seq_list and $rnai_seqlist created.\n";
print "\n";

# #################### get Interpro-derived gene associations for GO ####################

open (WORM2GO, "$worm2go_infile") || die "Wormpep2go file $worm2go_infile not found. $!\n";
open (GENE_ASSOC, ">$gene_assoc_outfile") || die "Couldn't open file $gene_assoc_outfile. $!\n";

open (GO_TERM_TO_ONTOLOGY_LIST, "$onto_hash_file") || die "Can't open $onto_hash_file -- $!\n";
open (CDS_TO_LOCUS_LIST, "$cds_to_locus_file") ||  die "Can't open $cds_to_locus_file -- $!\n";

# 1. In wormpep2go:
#
#     First column is CDS.  Save.
#         Use as gene name if locus not assigned.
#     Second column is WP number -- save, use as database ID number.
#     Third column is Interpro entry: keep (to modify evidence code).
#     Fourth column is GO number: keep.
#
# 2. Spit out, for each line of wormpep2go:
#
#     Column 1: WB  [Invariant database abbreviation.]
#     Column 2: WP:---  [ID number, from col. 2 of wormpep2go.]
#     Column 3: Output of CDS-to-gene hash, given CDS.  [Gene.]
#     Column 4: [NOT] column!    <-- must be sure this is in.
#     Column 5: GO term, taken straight.
#     Column 6: PUBMED:11159333  [Most recent Interpro reference.]
#     Column 7: IEA  [Evidence -- Inferred by Electronic Annotation.]
#     Column 8: [With] -- used for various things like IGI, IPI. <-- mandatory skip.
#     Column 9: Ontology letter (P, F, C): from GO hash lookup.
#     Column 10: CDS name, and gene name if it's distinct, separated by stave.
#                  Eventually we want classical synonyms in here too.

@go_term_to_ontology_list = <GO_TERM_TO_ONTOLOGY_LIST>;
@cds_to_locus_list = <CDS_TO_LOCUS_LIST>;

chomp(@go_term_to_ontology_list);
chomp(@cds_to_locus_list);

%go_term_to_ontology_hash = @go_term_to_ontology_list;
%cds_to_locus_hash = @cds_to_locus_list;

# Note: perl hashes appear to choke on strings containing "." or ":".

$whoami = `whoami`;
$uname = `uname -a`;
$date = `date`;
$hostname = `hostname`;

chomp($whoami);
chomp($uname);
chomp($date);
chomp($hostname);

$bugtracker_file = "bugtracking_on_".$orig_datestamp;
open (BUGTRACKER, ">$bugtracker_file");

print GENE_ASSOC "!version: \$Revision:  \$\n";
print GENE_ASSOC "!date: \$Date:  \$\n";
print GENE_ASSOC "! \n";
print GENE_ASSOC "! Gene assignments for Caenorhabditis elegans, Wormbase.\n";
print GENE_ASSOC "! IEA annotations were generated from ACeDB by Tablemaker.\n";
print GENE_ASSOC "! Processed for GO by $author_address.\n";
print GENE_ASSOC "! Processed for GO on $whoami\@$hostname account.\n";
print GENE_ASSOC "! Machine specs: $uname \n";
print GENE_ASSOC "! WB_Date: $date \n";
print GENE_ASSOC "! \n";

while (<WORM2GO>) 
{
    if ($_ =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+)/) 
    {
        $cds_literal = $1;
        $wp_id = $2;
        $interpro_code = $3;
        $go_term_literal = $4;
        $cds_hash_friendly = $cds_literal;
        $cds_hash_friendly =~ s/\./\_/g;
        $go_term_hash_friendly = $go_term_literal;
        $go_term_hash_friendly =~ s/\:/\_/g;

        # Print "DB" (1st) column.
        print GENE_ASSOC "WB";
        print GENE_ASSOC "\t";

        # Print "DB_Object_ID" (2cd) column -- i.e., "Gene_ID".
        print GENE_ASSOC "$wp_id";
        print GENE_ASSOC "\t";

        # Print "DB_Object_Symbol" (3rd) column -- i.e., "Gene_Symbol".
        print GENE_ASSOC "$cds_to_locus_hash{$cds_hash_friendly}";
        print GENE_ASSOC "\t";

        # Print "NOT" (4th) column.
        # In most cases this is pro forma, thus blank.
        print GENE_ASSOC "\t";

        # Print "GOid" (5th) column.
        print GENE_ASSOC "$go_term_literal";
        print GENE_ASSOC "\t";

        # Print "DB:Ref" (6th) column.
        print GENE_ASSOC "PUBMED:11159333";
        print GENE_ASSOC "\t";

        # Print "Evidence" (7th) column.
        print GENE_ASSOC "IEA";
        print GENE_ASSOC "\t";

        # Print "With" (8th) column.
        # In most cases this is pro forma, thus blank.  However, 
        #    it's a good place to put the Interpro code.
        print GENE_ASSOC "$interpro_code";
        print GENE_ASSOC "\t";

        # Print "Aspect" (9th) column.
        print GENE_ASSOC "$go_term_to_ontology_hash{$go_term_hash_friendly}";
        print GENE_ASSOC "\t";

        # Print "Name(|Name)" (10th) column.
        # In most cases this is pro forma, thus blank.
        print GENE_ASSOC "\t";

        # Print "Synonym(|Synonym)" (11th) column.
        if ($cds_literal eq $cds_to_locus_hash{$cds_hash_friendly}) 
        {
            print GENE_ASSOC "\t";
        }
        else 
        {
            print GENE_ASSOC "$cds_literal";
            print GENE_ASSOC "\t";
        }
        # Print "DB_Object_Type" (12th) column -- usually "protein".
        # But note that, for RNA products of genes, one needs "transcript".
        print GENE_ASSOC "protein";
        print GENE_ASSOC "\t";

        # Print "Taxon_ID" (13th) column -- always "taxon:6239".
        print GENE_ASSOC "taxon:6239";
        print GENE_ASSOC "\t";

        # Print "Date" (14th) column -- e.g., "20010127" for Jan. 27, 2002.
        print GENE_ASSOC "$datestamp";
        print GENE_ASSOC "\n";
    }
}

close WORM2GO;
close GENE_ASSOC;
close GO_TERM_TO_ONTOLOGY_LIST;
close CDS_TO_LOCUS_LIST;

# #################### add RNAi-derived gene associations for GO ####################

# Not needed?
# open (RNA2GO, "$rnai2go_file") || die "RNAi-to-GO file $rnai2go_file not opened or found. $!\n";

# Format of $rnai2go_file is:
#
# "JA:B0261.4"	"Bmd"	"B0261.4"	"WP:CE07706"
# "JA:B0261.4"	"Emb"	"B0261.4"	"WP:CE07706"
# "JA:B0261.4"	"Gro"	"B0261.4"	"WP:CE07706"

open (RNAI_SEQLIST, "$rnai_seqlist") 
    || die "RNAi-specific sequence list $rnai_seqlist could not be opened. $!\n";
open (GENE_ASSOC, ">>$gene_assoc_outfile") || die "Couldn't open file $gene_assoc_outfile. $!\n";

open (GO_TERM_TO_ONTOLOGY_LIST, "$onto_hash_file") || die "Can't open $onto_hash_file -- $!\n";
open (CDS_TO_LOCUS_LIST, "$cds_to_locus_file") ||  die "Can't open $cds_to_locus_file -- $!\n";
@go_term_to_ontology_list = <GO_TERM_TO_ONTOLOGY_LIST>;
@cds_to_locus_list = <CDS_TO_LOCUS_LIST>;
chomp(@go_term_to_ontology_list);
chomp(@cds_to_locus_list);
%go_term_to_ontology_hash = @go_term_to_ontology_list;
%cds_to_locus_hash = @cds_to_locus_list;

@sequence=<RNAI_SEQLIST>;
{
    foreach $sequence (@sequence)
    {
        chomp($sequence);
        $temp_rnai2go_file = "temp_" . $rnai2go_file . "_for_" . $sequence;
        system "cat $rnai2go_file | grep $sequence > $temp_rnai2go_file";
        open (TEMP_RNA2GO, "$temp_rnai2go_file") || die "Couldn't open file $temp_rnai2go_file. $!\n";
        while (<TEMP_RNA2GO>) 
        {

# Relevant format is, e.g.:
# "SA:yk345c8"	"WT"	"6R55.1"

            $inputline = $_;
            $inputline =~ s/\"//g;
            if ($inputline =~ /^(\S+\:\S+)\s+WT/)
            {
                print BUGTRACKER "Ignoring $1, a WT RNAi.\n";
            }

# Another relevant format is, e.g.:
# "KK:AH6.5"	"Emb"	"AH6.5"	"WP:CE26846"
# Problems I've seen:
# Problem input line: JA:F21C3.5	Emb	F21C3.5	WP:CE05680 
# Problem input line: JA:F21C3.5	Gro	F21C3.5	WP:CE05680 
# Problem input line: JA:F21C3.5	Unc	F21C3.5	WP:CE05680 
# Problem input line: JA:F21C3.6	WT	F21C3.6	WP:CE05681 
# Ignoring SA:yk261a7, a WT RNAi.
# Problem input line: JA:F21F12.1	WT	F21F12.1	WP:CE09546 
# Problem input line: JA:F21F3.1	WT	F21F3.1	WP:CE09531 


            elsif ($inputline =~ /^(\S+\:\S+)\s+(\S{3,})\s+(\w+\.\w+)(.+)/)
            {
                $rnai_name = $1;
                $phenotype = $2;
                $cds_literal = $3;
                $wp_id = $4;
                if ($wp_id =~ /.+(WP\:\w+).*/) 
                {
                    $wp_id = $1;
                }
                else 
                {
                    $wp_id = "";
                }
                print BUGTRACKER "Trying to handle RNAi $rnai_name, phenotype $phenotype, and CDS $cds_literal (ID no. \"$wp_id\")\n";
                $cds_literal =~ tr/a-z/A-Z/;
                $cds_hash_friendly = $cds_literal;
                $cds_hash_friendly =~ s/\./\_/g;
                $temp_pheno_go_list = "temp_" . $pheno_go_list . "_for_" . $phenotype;

                system "cat $pheno_go_list | grep $phenotype > $temp_pheno_go_list";
                open (TEMP_PHENO2GO, "$temp_pheno_go_list") || die "Couldn't open file $temp_pheno_go_list. $!\n";
                while (<TEMP_PHENO2GO>) 
                {
                    $phenotype_description_line = $_;
                    chomp($phenotype_description_line);
                    print BUGTRACKER "Managed to read $phenotype_description_line.\n";
                    @phenotype_term = split /\t/, $phenotype_description_line;
                    foreach $phenotype_term (@phenotype_term) 
                    {
                        if ($phenotype_term =~ /.*(GO\S+)/) 
                        {
                            $phenotype_term = $1;
                            $phenotype_term_hash_friendly = $phenotype_term;
                            $phenotype_term_hash_friendly =~ s/\:/\_/g;
        
                            # Print "DB" (1st) column.
                            print GENE_ASSOC "WB";
                            print GENE_ASSOC "\t";
        
                            # Print "DB_Object_ID" (2cd) column -- i.e., "Gene_ID".
                            print GENE_ASSOC "$wp_id";
                            print GENE_ASSOC "\t";
        
                            # Print "DB_Object_Symbol" (3rd) column -- i.e., "Gene_Symbol".
                            print GENE_ASSOC "$cds_to_locus_hash{$cds_hash_friendly}";
                            print GENE_ASSOC "\t";
        
                            # Print "NOT" (4th) column.
                            # In most cases this is pro forma, thus blank.
                            print GENE_ASSOC "\t";

                            # Print "GOid" (5th) column.
                            print GENE_ASSOC "$phenotype_term";
                            print GENE_ASSOC "\t";
            
                            # Print "DB:Ref" (6th) column.
                            # If $rnai_name starts with JA: -- PMID:11099033
                            # eventually, for 147 of these: -- PMID:11483502
                            # If $rnai_name starts with KK: -- PMID:11137018
                            # If $rnai_name starts with SA: -- PMID:11231151
                            # If $rnai_name starts with TH: -- PMID:11099034
                            # If $rnai_name starts with [cgc*]: -- WB:cgc*
                            # E.g., [cgc4722]: -- WB:cgc4722
                            # If $rnai_name starts with [med*]: -- WB:med*
                            # E.g., [med8322073]: -- WB:med8322073

                            if ($rnai_name =~ /JA\:.+/) 
                            {
                                print GENE_ASSOC "PMID\:11099033";
                                print GENE_ASSOC "\t";
                            }
                            elsif ($rnai_name =~ /KK\:.+/)
                            {
                                print GENE_ASSOC "PMID\:11137018";
                                print GENE_ASSOC "\t";
                            }
                            elsif ($rnai_name =~ /SA\:.+/)
                            {
                                print GENE_ASSOC "PMID\:11231151";
                                print GENE_ASSOC "\t";
                            }
                            elsif ($rnai_name =~ /TH\:.+/)
                            {
                                print GENE_ASSOC "PMID\:11099034";
                                print GENE_ASSOC "\t";
                            }
                            elsif ($rnai_name =~ /^\[(cgc\d+)\]\:.+/)
                            {
                                print GENE_ASSOC "WB\:$1";
                                print GENE_ASSOC "\t";
                            }
                            elsif ($rnai_name =~ /^\[(med\d+)\]\:.+/)
                            {
                                print GENE_ASSOC "WB\:$1";
                                print GENE_ASSOC "\t";        
                            }
                            else
                            {
                                print GENE_ASSOC "WB\:$rnai_name";
                                print GENE_ASSOC "\t";
                            }
        
                            # Print "Evidence" (7th) column.
                            print GENE_ASSOC "IMP";
                            print GENE_ASSOC "\t";
        
                            # Print "With" (8th) column.
                            # In most cases this is pro forma, thus blank.  However,
                            #    it's a good place to put the RNAi.
                            print GENE_ASSOC "WB\:$rnai_name";
                            print GENE_ASSOC "\t";

                            # Print "Aspect" (9th) column.
                            print GENE_ASSOC "$go_term_to_ontology_hash{$phenotype_term_hash_friendly}";
                            print GENE_ASSOC "\t";
        
                            # Print "Name(|Name)" (10th) column.
                            # In most cases this is pro forma, thus blank.
                            print GENE_ASSOC "\t";
        
                            # Print "Synonym(|Synonym)" (11th) column.
                            if ($cds_literal eq $cds_to_locus_hash{$cds_hash_friendly})
                            {
                                print GENE_ASSOC "\t";
                            }
                            else
                            {
                                print GENE_ASSOC "$cds_literal";
                                print GENE_ASSOC "\t";
                            }

                            # Print "DB_Object_Type" (12th) column -- always "protein".
                            print GENE_ASSOC "protein";
                            print GENE_ASSOC "\t";
        
                            # Print "Taxon_ID" (13th) column -- always "taxon:6239".
                            print GENE_ASSOC "taxon:6239";
                            print GENE_ASSOC "\t";

                            # Print "Date" (14th) column -- e.g., "20010127" for Jan. 27, 2002.
                            print GENE_ASSOC "$datestamp";
                            print GENE_ASSOC "\n";

                        } # End "if ($phenotype_term =~ /^GO\w+/)".
                    } # End "foreach $phenotype_term (@phenotype_term)".
                } # End "while (<TEMP_PHENO2GO>)".
                system "rm $temp_pheno_go_list";
            }
            else 
            {
                chomp($inputline);
                print BUGTRACKER "Problem input line: $inputline \n";
            }
        }
        system "rm $temp_rnai2go_file";
    }
}

close TEMP_RNA2GO;
close RNAI_SEQLIST;
close GENE_ASSOC;
close GO_TERM_TO_ONTOLOGY_LIST;
close CDS_TO_LOCUS_LIST;

# #################### convert entire gene association file into an .ace file ####################

open (SEQLIST, "$seq_list") || die "Sequence list file $seq_list could not be opened. $!\n";
open (GO2SEQ_ACE, ">$go2seq_ace") || die "Couldn't open file $go2seq_ace. $!\n";

print GO2SEQ_ACE "\/\/ These data were generated by $author_address\n";
print GO2SEQ_ACE "\/\/ on $date, using make_GO_Web_and_ace_file.pl.\n";
print GO2SEQ_ACE "\/\/ The machine specs were:\n";
print GO2SEQ_ACE "\/\/ $uname\n";
print GO2SEQ_ACE "\n";

# Format of GENE_ASSOC will be:
#
# WB	WP:CE20469	B0564.1		GO:0000175	PUBMED:11159333	IEA	INTERPRO:IPR001247	F			protein	taxon:6239
# WB	WP:CE20469	B0564.1		GO:0003723	PUBMED:11159333	IEA	INTERPRO:IPR001247	F			protein	taxon:6239
# [...]
# WB	WP:CE20469	B0564.1		GO:0002119	SA:yk198e11	IMP	SA:yk198e11	P			protein taxon:6239

# Format of GO2SEQ_ACE will be:
#
# // data dumped from tree display
#
# Sequence : "B0564.1"
# GO_term	 "GO:0000175"	\\ IEA: from INTERPRO:IPR001247
# GO_term	 "GO:0003723"	\\ IEA: from INTERPRO:IPR001247
# [...]
# GO_term	 "GO:0002119"	\\ IMP: from SA:yk198e11

# printed out with something like:
#        print GO2SEQ_ACE "Sequence : \"$cds_literal\"\n";
#        print GO2SEQ_ACE "GO_term\t \"$go_term_literal\"   \/\/ IEA: from $interpro_code\n";
#        print GO2SEQ_ACE "\n";


# was:
# open (GENE_ASSOC, "$gene_assoc_outfile") || die "Couldn't open file $gene_assoc_outfile. $!\n";

open (ROUGH_PROBLEM_NAMES, ">>$problem_name_file.1") || die "Couldn't open rough file $problem_name_file.1. $!\n";
open (SHODDY_INPUT_LINES, ">>$shoddy_input_file") || die "Couldn't open shoddy input lines file $shoddy_input_file. $!\n";

@sequence=<SEQLIST>;
{ 
    foreach $sequence (@sequence) 
    {
        chomp($sequence);
        print GO2SEQ_ACE "\n";
        print GO2SEQ_ACE "Sequence : \"$sequence\"\n";
        $temp_assoc_file = "temp_" . $gene_assoc_outfile . "_for_" . $sequence;
        system "cat $gene_assoc_outfile | grep $sequence > $temp_assoc_file";
        open (TEMP_GENE_ASSOC, "$temp_assoc_file") || die "Couldn't open file $gene_assoc_outfile. $!\n";

# Needs to be able to handle two different inputs:
# WB	WP:CE01453	sra-8		GO:0007165	PUBMED:11159333	IEA	INTERPRO:IPR000344	P		AH6.12	protein	taxon:6239
# WB	WP:CE07654	AH6.13		GO:0000155	PUBMED:11159333	IEA	INTERPRO:IPR000344	F			protein	taxon:6239

# Have been choking on:
# WB	WP:CE01450	gcy-1		GO:0004672	PUBMED:11159333	IEA	INTERPRO:IPR000719	F		AH6.1	protein	taxon:6239
# and...
# WB	WP:CE21025	air-1		GO:0007345	WB:SA:yk364b4	IMP	WB:SA:yk364b4	P		K07C11.2	protein	taxon:6239
# WB	WP:CE21025	air-1		GO:0007276	WB:SA:yk364b4	IMP	WB:SA:yk364b4	P		K07C11.2	protein	taxon:6239
# WB	WP:CE21025	air-1		GO:0007320	WB:SA:yk364b4	IMP	WB:SA:yk364b4	P		K07C11.2	protein	taxon:6239
# WB	WP:CE21025	air-1		GO:0007338	WB:SA:yk364b4	IMP	WB:SA:yk364b4	P		K07C11.2	protein	taxon:6239

        while (<TEMP_GENE_ASSOC>) 
        {
            $input_line = $_;
# Was:
#            if ($input_line =~ /^WB\s+WP:CE\d+\W+(\S+\.\S+)\W+(GO:\d+)\W+\w+:\w+\W+(\w{3})\W+(\w+:\w+)\W+.+/) 
            if ($input_line =~ /^WB\s+WP:CE\d+\W+(\S+\.\S+).+(GO:\d+)\W+\S+\W+(\w{3})\W+(\S+)\W+.+/)
            {
                $cds_literal=$1;
                $go_term_literal=$2;
                $evidence_code=$3;
                $evidence_reference=$4;
                if ($sequence eq $cds_literal) 
                {
                    print GO2SEQ_ACE "GO_term\t \"$go_term_literal\"   \/\/ $evidence_code: from $evidence_reference\n";
                }
            }
            elsif ($input_line =~ /^WB\s+WP:CE\d+\W+\S+\W+(GO:\d+)\W+\S+\W+(\w{3})\W+(\S+)\W+\w{1}\W+(\S+\.\S+)/)
            {
                $go_term_literal=$1;
                $evidence_code=$2;
                $evidence_reference=$3;
                $cds_literal=$4;
                if ($sequence eq $cds_literal)
                {
                    print GO2SEQ_ACE "GO_term\t \"$go_term_literal\"   \/\/ $evidence_code: from $evidence_reference\n";
                }
            }
            elsif ($input_line =~ /^WB\s+WP:CE\d+\W+(\S+)\W+(GO:\d+)/) 
            {
                print ROUGH_PROBLEM_NAMES "$1\n";
            }
            else 
            {
                print SHODDY_INPUT_LINES "$input_line";
            }
        }
        close TEMP_GENE_ASSOC;
        system "rm $temp_assoc_file";
    }
}

close ROUGH_PROBLEM_NAMES;
close SHODDY_INPUT_LINES;
system "sort $problem_name_file.1 | uniq - > $problem_name_file";
system "rm $problem_name_file.1";
