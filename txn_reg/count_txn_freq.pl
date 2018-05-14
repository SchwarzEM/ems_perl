#!/usr/bin/perl

# count_txn_freq.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/29/2007.
# Purpose: just list transcription factors in table by most frequently named, downward.
# Example: ./count_txn_freq.pl txn_fac1_v20.txt > txn_fac_freqs_29oct2007.txt

use strict;
use warnings;

my %names = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) /xms ) { 
        $names{$1}++;
    }
}

my @namelist = sort { $names{$b} <=> $names{$a} } sort keys %names;

foreach my $name (@namelist) { 
    print "$name\t$names{$name}\n";
}

