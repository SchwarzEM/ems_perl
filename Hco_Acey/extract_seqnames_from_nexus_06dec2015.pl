#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    while ( $input =~ /\' ([^\'\/]+) \/ /gxms ) { 
        print "$1\n";
    }
}


