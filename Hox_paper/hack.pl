#!/usr/bin/perl

use strict;
use warnings;

while (my $input = <>) { 
    if ( ($input =~ / \D (\d+) \s+ aa /xms ) and ($1 >= 100) ) { 
        print $input;
    }
}
