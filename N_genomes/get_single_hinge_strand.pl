#!/usr/bin/env perl
# get_single_hinge_strand.pl -- Erich Schwarz <ems394@cornell.edu>, 3/8/2018.

use strict;
use warnings;
use autodie;

my $print = 0;;

while (my $input = <>) {
    chomp $input;

    # make decision whether to be in printing mode or not
    if ( $input =~ /\A [>] \S+ (\d+) \z/xms ) { 
        my $index = $1;
        $index    = ($index / 2);
        if ( $index == int($index) ) {
            $print = 1;
        }
        else {
            $print = 0;
        }
    }

    if ($print) {
        print "$input\n";
    }
}
