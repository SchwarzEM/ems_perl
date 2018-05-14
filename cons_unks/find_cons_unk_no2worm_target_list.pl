#!/usr/bin/env perl

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input !~ / .* c_elegans\| .+ c_elegans\| /xms ) { 
        print "$input\n";
    }
}


