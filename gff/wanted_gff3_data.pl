#!/usr/bin/env perl

# wanted_gff3_data.pl -- Erich Schwarz <emsch@its.caltech.edu>, 5/29/2009.
# Purpose: filter out GFF3s, or, see if a mapped GFF3 has 'extra' lines.

use strict;
use warnings;

while (my $input = <>) { 
    if (     ( $input =~ /\A \S+ \t Coding_transcript \t (coding_region_of_exon|CDS|gene|mRNA) \t /xms ) 
          or ( $input =~ /\A [#] /xms ) ) { 
        print $input;
    }
}

