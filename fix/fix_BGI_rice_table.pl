#!/usr/bin/perl

# fix_BGI_rice_table.pl

# Purpose: reformat Beijing GI rice table from tab-delimited ASCII to FASTA.

# Input format is:
#
# (1) . Functional classification
# (2).  Gene name in BGI database
# (3).  FGeneSH protein size [a.a.]
# (4).  Extent of hit (TblastN)
# (5).  AA identity (TblastN)
#
# (6).  Hits per gene (TblastN)
# (7).  Extent of hit (BlastP)
# (8).  AA identity (BlastP)
# (9).  Rice accession number
# (10). Identity with rice cDNA
#
# (11). Arab accession number
# (12). Identity with arab cDNA
# (13). Gene description (from rice or arab)
# (14). Protein sequence for rice prediction
#
# The format used for FASTA headers was as follows
# (all of these are on the header line, with three spaces between each):
#
# >Rice_BGI_[Gene name in BGI database]                == (2) --> array '1', since first no. is '0'.
# Size: [FGeneSH protein size [a.a.]]                  == (3) --> array '2'
# Rice acc. no.: [Rice accession number or "none"]     == (9) --> array '8'
# Arab. acc. no.: [Arab accession number or "none"]    == (11) --> array '10'
# Oryza sativa L. ssp. indica, cultivar 93-11; PMID:11935017.
#
# The protein sequences were reformatted to have <= 60 residues per line 
#     -- protein sequence == (14) -> array '13'

use strict;

print "Input file?   ";
chomp (my $input_file = <STDIN>);
my $output_file = $input_file . ".fasta";

open (INPUT, "$input_file")    || die;
open (OUTPUT, ">$output_file") || die;

while (<INPUT>) 
{
    chomp (my $input_line = $_);
    my @input_array = split /\t/, $input_line;
    print OUTPUT ">Rice_BGI_";
    print OUTPUT "$input_array[1]";
    print OUTPUT "   ";
    print OUTPUT "$input_array[2] residues";
    print OUTPUT "   ";
    unless ($input_array[8] =~ /^\s*$/) { print OUTPUT "[Rice acc.: $input_array[8]; " };
    if ($input_array[8] =~ /^\s*$/) { print OUTPUT "[Rice acc.: none; " };
    unless ($input_array[10] =~ /^\s*$/) { print OUTPUT "Arab. acc.: $input_array[10].]  " };
    if ($input_array[10] =~ /^\s*$/) { print OUTPUT "Arab. acc.: none.]  " };
    print OUTPUT "O. sativa ssp. indica (93-11); PMID:11935017\n";
    my $line_to_split = $input_array[13];
    $line_to_split =~ s/\s//g;
    while ($line_to_split =~ /^([\w]{60})(\w*)/) 
    {
        my $line_to_print = $1;
        print OUTPUT "$line_to_print\n";
        $line_to_split = $2;
    }
    print OUTPUT "$line_to_split\n";
}

close INPUT;
close OUTPUT;
