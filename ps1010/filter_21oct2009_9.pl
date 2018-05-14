#!/usr/bin/env perl

# filter_21oct2009_9.pl -- get various alternative prot-genes (from pre-filtered "Transcript \"[^\"\s]+\" lines from huge GFF files).

use strict;
use warnings;

my %OK_terms = ( 'GeneMarkHMM' => 1,
                 'Genefinder' => 1,
                 'jigsaw' => 1,
                 'mSplicer_orf' => 1,
                 'mSplicer_transcript' => 1,
                 'twinscan' => 1, );

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ / \A \S+ \t (\S+) /xms ) { 
        my $term = $1;
        if ( $OK_terms{$term} ) { 
            print "$input\n";
        }
    }
} 

