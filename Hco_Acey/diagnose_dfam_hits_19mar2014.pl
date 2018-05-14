#!/usr/bin/env perl

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input !~ /\A [#] /xms ) { 
        my @fields = split /\s+/, $input;
        if( $fields[12] <= 1e-20 ) { 
            print "$fields[0]\t$fields[1]\t$fields[12]\n";
        }
    }
}

