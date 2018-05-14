#!/usr/bin/env perl

# censor_UCSC_worm_wigs.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/11/2010.
# Purpose: given a .wig that may have spikes in it, censor everything that's not a comment line or chr[X] sequence.

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if (     ( $input =~ /\A chr(?:I|II|III|IV|V|X|M) /xms ) 
          or ( $input =~ /\A track [ ] type=bedGraph [ ] name=\"[^\"\s]+\" /xms ) ) { 
        print "$input\n";
    }
}

