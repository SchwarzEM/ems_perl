#!/usr/bin/env perl

use strict;
use warnings;

while (my $input = <>) {
    chomp $input;
    # Acey_2012.08.05_0008.g153.t3  -          2OG-FeII_Oxy         PF03171.15   2.1e-07   33.2 
    if ( ( $input !~ /\A [#] /xms ) and ( $input =~ /\A (?: \S+ \s+){3} (\S+) \s+ (\S+) /xms ) ) { 
        my $pfam_acc = $1;
        my $e_value  = $2;
        if ( $e_value <= 1e-06 ) { 
            print "$pfam_acc\n";
        }
    }
}


