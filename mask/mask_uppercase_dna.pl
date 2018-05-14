#!/usr/bin/perl -w

# mask_uppercase_dna.pl -- Erich Schwarz, 11/4/07, <emsch@its.caltech.edu>.
# Purpose: use either 'X' or 'N' to mask uppercase nucleotides in a FASTA file.

use strict;
use warnings;

if ( $#ARGV != 1) { 
    die "Format: ./mask_uppercase_dna.pl [X|N (as mask nt)] [input ASCII text file]\n";
}

my $mask = shift @ARGV;
if ( ($mask ne 'X') and ($mask ne 'N') ) {
    die "Mask nucleotide incorrectly specified as: $mask\n";
}

while (my $input = <>) {
    chomp $input;
    unless ($input =~ /^>/) {
        $input =~ s/[A-Z]/$mask/g;
    }
    print "$input\n";
}

