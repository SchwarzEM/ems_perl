#!/usr/bin/perl -w

# clean_up_YeastPep.pl
# Erich Schwarz <emsch@its.caltech.edu>, 11/27/02.

# Purpose: Strip "ORFP:" out of ">ORFP:[etc.]" FASTA headers of YeastPep.
# 
# YeastPep == orf_trans.fasta.gz  1908 KB  07/25/2002  05:28:00 PM
# from        ftp://genome-ftp.stanford.edu/pub/yeast/data_download/sequence/genomic_sequence/orf_protein/

print "input file: ";
$input = <STDIN>;
chomp($input);
$output = "YeastPep__which_was__" . $input;

open (INPUT, "$input") || die;
open (OUTPUT, ">$output") || die;

while (<INPUT>) 
{
    chomp ($input_line = $_);
    if ($input_line =~ />ORFP:(\S+)(.*)/) 
    {
        print OUTPUT ">$1    $2\n";
    }
    elsif ($input_line =~ /(\S+)\*$/)
    {
        print OUTPUT "$1\n";    # finally fixed bug here, 11/27/02!
    }
    else 
    {
        print OUTPUT "$input_line\n";
    }
}

close INPUT;
close OUTPUT;
