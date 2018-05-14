#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $source = q{};
my $prefix = q{};

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A ==> [ ] proteomes \/ ([A-Za-z_]+) \. /xms ) { 
        $source = $1;
    }
    elsif ( $input =~ /\A > ([A-Za-z]+) _ /xms ) { 
        $prefix = $1;
        print "    '$prefix' => '$source',\n";
    }
}

