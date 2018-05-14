#!/usr/bin/env perl

use strict;
use warnings;

while (my $input = <>) { 
    my $tabcount = ( $input =~ tr/\t/\t/ );
    print "$tabcount\n";
}


