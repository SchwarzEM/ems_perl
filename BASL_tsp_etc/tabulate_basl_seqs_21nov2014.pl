#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    my $species = 'TBD';
    if ( $input =~ /\A > (\S+) \s+ (.+) \z/xms ) { 
        my $name = $1;
        my $data = $2;
        if ( $data =~ /\A (.+) \s+ \[ (.+) \] \s* /xms ) { 
            $species = $2;
            $data    = $1;
        }
        elsif ( $name =~ /\A Aquca_ /xms ) { 
            $species = 'Aquilegia coerulea';
        }
        elsif ( $name =~ /\A Spipo /xms ) {
            $species = 'Spirodela polyrhiza';
        } 
        print "$name\t$data\t$species\n";
    }
}

