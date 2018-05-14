#!/usr/bin/perl

# DIANA-1.1.2.pl

# by Melissa Michelitsch,
# with modifications by Erich Schwarz <emsch@its.caltech.edu>, 1/18/02.
#
# Derived from the original DIANA, by:
#   Melissa Michelitsch <mmichel@itsa.ucsf.edu>
#       [e-mail address on 1/18/02]
#   Melissa Michelitsch <melissa_michelitsch@chiron.com>
#       [e-mail address after UCSF address stops being used]
# 
# DIANA was described in the article:
# 
# Michelitsch and Weissman (2000).  A census of 
# glutamine/asparagine-rich regions: implications for
# their conserved function and the prediction of novel prions.
# Proc. Natl. Acad. Sci. U.S.A. vol. 97, pp. 11910-11915.
# 
# DIANA-1.1.pl:   1/18/02.  Changes to make DIANA easier to run.
# DIANA-1.1.1.pl: 1/30/02.  Bugfix, to correctly parse FASTA files.
# DIANA-1.1.2.pl: 1/30/02.  Change to make hitlists simpler/more useful.
# 

print "What protein FASTA file is to be scanned? ";
$protein_fasta_input = <STDIN>;
chomp ($protein_fasta_input);

$diana_output           = $protein_fasta_input . ".diana-output";
$filtered_diana_output  = $protein_fasta_input . ".diana-filtered-output";
$diana_seq              = $protein_fasta_input . ".diana-seq";
$diana_hitlist          = $protein_fasta_input . ".diana-hitlist";
$diana_hashable_hitlist = $protein_fasta_input . ".diana-hashable-hitlist";

open (DIANA_OUTPUT, ">$diana_output") || die "Can't open DIANA output file $diana_output. $! \n";

# Generate reformatted versions of a given protein FASTA database.

open (DIANA_SEQ, ">$diana_seq") || die ("Can't open DIANA-formatted sequence file $diana_seq. $! \n"); 

$reformatting_underway = "no";

open (PROTEIN_FASTA_INPUT, "$protein_fasta_input") || die ("Can't open input file $protein_fasta_input. $!\n");

while (<PROTEIN_FASTA_INPUT>) 
{
    $input_line = $_;
    chomp ($input_line);
    if ($input_line =~ /^>gi\|(\w+)\|/ && $reformatting_underway eq "no") 
    {
        print DIANA_SEQ ">$1\n";
        $reformatting_underway = "yes";
    }
    elsif ($input_line =~ /^>(\S+)\s+/ && $reformatting_underway eq "no") 
    {
        print DIANA_SEQ ">$1\n";
        $reformatting_underway = "yes";
    }
    elsif ($input_line =~ /^>(\S+)$/ && $reformatting_underway eq "no")
    {
        print DIANA_SEQ ">$1\n";
        $reformatting_underway = "yes";
    }


    elsif ($input_line =~ /^>gi\|(\w+)\|/ && $reformatting_underway eq "yes") 
    {
        print DIANA_SEQ "\*";   
        print DIANA_SEQ "\n";   
        print DIANA_SEQ ">$1\n";
    } 
    elsif ($input_line =~ /^>(\S+)\s+/ && $reformatting_underway eq "yes")
    {
        print DIANA_SEQ "\*";
        print DIANA_SEQ "\n";
        print DIANA_SEQ ">$1\n";
    }
    elsif ($input_line =~ /^>(\S+)$/ && $reformatting_underway eq "yes")
    {
        print DIANA_SEQ "\*";
        print DIANA_SEQ "\n";
        print DIANA_SEQ ">$1\n";
    }
    else 
    {
        print DIANA_SEQ "$input_line";
    }
}
close PROTEIN_FASTA_INPUT;
close DIANA_SEQ;

# Now, proceed to do the analysis itself.

unless (open (PROTEIN_FASTA_INPUT, "$diana_seq"))
{
    die ("Can't open DIANA-formatted protein file $diana_seq. $!\n");
}

## Make title array ## 

$x=0; 
$seqnum=1; 

while ($line = <PROTEIN_FASTA_INPUT>) 
{ 
    chop $line; 
    if ($line =~ /^>/) 
    { 
        $title[$x] = $line;  
        $x++;
    }
} 

$titlecnt = @title;

## Make sequence array ##

unless (open (EDIT, "$diana_seq")) 
{
    die ("Won't open.\n"); 
}
 
$_ = join '',<EDIT>;     # Concatenate all input lines #
s/>.*//gm;               # Delete title lines (e.g. ">Contig2")
s/\s//g;                   # Delete white space 

@array = (split(/\*/,$_)); 

##Make Letter-sequence array ## 

$b=0; 

while ($title[$b] ne "")
{
    $_ = $array[$b];
    s/\s//gm; 
    @letarray = (split(//,$_)); 
    $ltrcount = 0;
    $ltrcount += @letarray; 

    $maxQ=0;

    if ($ltrcount < 80) 
    {
        $x=($ltrcount-1);
    } 
    else 
    {
        $x=79;
    } 
    $y=0; 
    $z=0; 
    $r=0; 


##make 80-mers## 

    $ATmer = join("", @letarray[$z..$x]); 

##split seq##

    @mer = (split(/Q|N|q|n/,$ATmer));

    if ($letarray[$x] !~ /Q|N|q|n/) 
    {
        $numQ = (@mer-1); 
    }
    else 
    {
        $numQ = (@mer); 
    }

## ref. seq. by # of Q's ## 

    while ($letarray[$x] ne "") 
    {
        if ($z != 0) 
        { 
            if ($letarray[$x] =~ /Q|N|q|n/) 
            { 
                $numQ++;
            } 
            $k = ($z-1);

            if ($letarray[$k] =~ /Q|N|q|n/) 
            {
                $numQ--;
            }
        }

##determine max Q ## 

        if ($maxQ <= $numQ) 
        { 
            $ATmer = join("", @letarray[$z..$x]); 
            $Qarray{$numQ} = $ATmer; 
            ($maxQ = $numQ) ##write to assoc. array by max number Q calling up assoc. seq. ## 
        }
        $x++; 
        $z++;
        $r++;
    }

print DIANA_OUTPUT "$title[$b]\t"; 
print DIANA_OUTPUT "$Qarray{$maxQ}\t"; 
$seqstr = $Qarray{$maxQ}; 
@Qcnt = (split(/Q|q/,$seqstr)); 
$numQ = (@Qcnt-1); 
@Ncnt = (split(/N|n/,$seqstr));
$numN = (@Ncnt-1);

##Output##

print DIANA_OUTPUT "$ltrcount\t$maxQ\t$numQ\t$numN\n";
$b++;

} 

close PROTEIN_FASTA_INPUT;
close EDIT;
close DIANA_OUTPUT;

open (DIANA_OUTPUT, "$diana_output") 
    || die "Can't open raw DIANA output file $diana_output! $! \n";
open (FILTERED_DIANA_OUTPUT, ">$filtered_diana_output") 
    || die "Can't open filtered DIANA output file $filtered_diana_output! $! \n";
open (DIANA_HITLIST, ">$diana_hitlist") 
    || die "Can't open DIANA hitlist file $diana_hitlist! $! \n";
open (DIANA_HASHABLE_HITLIST, ">$diana_hashable_hitlist") 
    || die "Can't open DIANA hashable hitlist file $diana_hashable_hitlist! $! \n";

while (<DIANA_OUTPUT>) 
{
    $input_line = $_;
    chomp ($input_line);
    foreach ($input_line) 
    {
        if ($input_line =~ /^>(.+)\s+.+\t\d+\t(\d+)\t\d+\t\d+/) 
        {
            $diana_hit_name = $1;
            $score_number = $2;
            if ($score_number > 29) 
            {
                print FILTERED_DIANA_OUTPUT "$input_line\n";
                print DIANA_HITLIST "$diana_hit_name\n";
                $diana_hit_name =~ s/\./_/g;
                print DIANA_HASHABLE_HITLIST "$diana_hit_name\n";
            }
        }
    }
}

close DIANA_OUTPUT;
close FILTERED_DIANA_OUTPUT;
close DIANA_HITLIST;
close DIANA_HASHABLE_HITLIST;
