#!/usr/bin/env perl

# wanted_gff2_data.pl -- Erich Schwarz <emsch@its.caltech.edu>, 5/28/2009.
# Purpose: pull out that subset of WormBase GFF2s which I currently think I know how to map to GFF3s.

use strict;
use warnings;

while (my $input = <>) { 
    if (     ( $input =~ /\A \S+ \t Coding_transcript \t (coding_exon|protein_coding_primary_transcript) \t /xms ) 
          or ( $input =~ /\A \S+ \t curated \t CDS \t /xms ) 
          or ( $input =~ /\A \S+ \t gene \t gene \t /xms )
          or ( $input =~ /\A [#] /xms ) ) { 
        print $input;
    }
}

