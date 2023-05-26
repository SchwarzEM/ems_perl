#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $name = q{};
my $type = q{};
my $desc = q{};

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\APosition_Matrix \s+ ["] (\S+) ["]/xms ) {
        $name = $1;
        $type = q{};
        $desc = q{};
    }
    elsif ( $input =~ /\AWeight/xms ) {
        $type = 'Weight';
    }
    elsif ( $input =~ /\ARemark \s+ (["] .+ ["])/xms ) {
       	$desc = $1;
        if ( $type eq 'Weight' ) {
            print "$name\t$type\t$desc\n";
        }
        $name = q{};
        $type = q{};
        $desc = q{};
    }
}

