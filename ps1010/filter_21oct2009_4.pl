#!/usr/bin/env perl

# filter_21oct2009_4.pl -- get an mRNA subset (from pre-filtered "Transcript \"[^\"\s]+\" lines from huge GFF files).

use strict;
use warnings;

my %OK_terms = ( 'Coding_transcript' => 1, );

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ / \A \S+ \t (\S+) /xms ) { 
        my $term = $1;
        if ( $OK_terms{$term} ) { 
            print "$input\n";
        }
    }
} 

