#!/usr/bin/env perl

# ucsc_ify_gffs.pl -- Erich Schwarz <emsch@its.caltech.edu>, 6/9/2009
# Purpose: put on dimwitted 'chr' prefix onto GFF sequence identifiers.

use strict;
use warnings;

while (my $input = <>) { 
    if ( ( $input !~ /\A \# /xms ) and ( $input =~ /\A \w /xms ) ) { 
        print 'chr';
    }
    print $input;
}


