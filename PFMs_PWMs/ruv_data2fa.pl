#!/usr/bin/env perl

# ruv_data2fa.pl -- Erich Schwarz <emsch@its.caltech.edu>, 4/8/2010.
# Purpose: convert Ruvinsky-provided data into a FASTA (for making into a PFM).

use strict;
use warnings;

my $i    = 0;
my $seq  = q{};
my %seqs = ();

while (my $input = <>) { 
    chomp $input;
    $input =~ s/\s*\z//;
    if ( $input =~ / \s ([ACGT]+) \z/xms ) { 
        $seq = $1;
        $i++;
        $seqs{$i} = $seq;
    }
}

foreach my $seq_no ( sort { $a <=> $b } keys %seqs ) { 
    print '>seq_', "$seq_no\n", ;
    print "$seqs{$seq_no}\n";
}

