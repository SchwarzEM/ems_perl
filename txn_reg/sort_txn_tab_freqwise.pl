#!/usr/bin/perl

# sort_txn_tab_freqwise.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/29/2007.
# Purpose: sort transcription factor table by most frequently named, downward.
# Example: ./sort_txn_tab_freqwise.pl txn_fac1_v20.txt > txn_fac1_byfreq_29oct2007.txt

use strict;
use warnings;

my %names             = ();
my %names2entries_ref = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A ( (\S+) .+) /xms ) { 
        $names{$2}++;
        push @{ $names2entries_ref{$2} }, $1;
    }
}

my @namelist = sort { $names{$b} <=> $names{$a} } sort keys %names;

foreach my $name (@namelist) { 
    foreach my $entry ( sort @{ $names2entries_ref{$name} } ) {
        print "$entry\n";
    }
}

