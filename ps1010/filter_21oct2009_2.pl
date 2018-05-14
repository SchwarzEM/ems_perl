#!/usr/bin/env perl

# filter_21oct2009_2.pl -- get a ncRNA subset (from pre-filtered "Transcript \"[^\"\s]+\" lines from huge GFF files).

use strict;
use warnings;

my %OK_terms = ( 'Non_coding_transcript' => 1,
                 'miRNA'                 => 1,
                 'ncRNA'                 => 1,
                 'rRNA'                  => 1,
                 'scRNA'                 => 1,
                 'snRNA'                 => 1,
                 'snlRNA'                => 1,
                 'snoRNA'                => 1,
                 'tRNA'                  => 1,         # added 10/27/2009, to cope with MtDNA
                 'tRNAscan-SE-1.23'      => 1, );

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ / \A \S+ \t (\S+) /xms ) { 
        my $term = $1;
        if ( $OK_terms{$term} ) { 
            print "$input\n";
        }
    }
} 

