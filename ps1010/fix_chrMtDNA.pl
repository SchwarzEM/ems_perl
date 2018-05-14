#!/usr/bin/env perl

# fix_chrMtDNA.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/11/2009.
# Purpose: change 'chrMtDNA' into 'chrM', which UCSC requires for worm GFFs.

use strict;
use warnings;

while (my $input = <>) { 
    if ( $input =~ /\A chrMtDNA \s /xms ) { 
        $input =~ s/chrMtDNA(\s.*)/chrM$1/;
    }
    print $input;
}

