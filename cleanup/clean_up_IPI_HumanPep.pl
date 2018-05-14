#!/usr/bin/perl -w

# clean_up_PfalPep.pl
# Erich Schwarz <emsch@its.caltech.edu>, 11/27/02.

# Purpose: Reformat IPI HumanPep headers so that simple protein names are clearly used.
# IPI HumanPep == ipi.HUMAN.fasta.gz 
# from       ftp://ftp.ebi.ac.uk/pub/databases/IPI/current/

# Sample header input:
#
# >IPI:IPI00047137.2|REFSEQ_XP...
#
# Sample desired header output:
# 
# >IPI00047137.2    REFSEQ_XP...

print "input file: ";
$input = <STDIN>;
chomp($input);
$output = $input . "_with_reformatted_headers";

open (INPUT, "$input") || die;
open (OUTPUT, ">$output") || die;

while (<INPUT>) 
{
    chomp ($input_line = $_);
    if ($input_line =~ />IPI:([^|]+)\|(.*)/) 
    {
        print OUTPUT ">$1    \|$2\n";
    }
    else 
    {
        print OUTPUT "$input_line\n";
    }
}

close INPUT;
close OUTPUT;
