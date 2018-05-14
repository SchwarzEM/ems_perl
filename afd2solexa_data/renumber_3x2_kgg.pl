#!/usr/bin/env perl

# renumber_3x2_kgg.pl -- Erich Schwarz <emsch@its.caltech.edu>, 7/24/2011.
# Purpose: get simple numbers from a hybrid cluster so I can map groups.

use strict;
use warnings;

my $front = q{};
my $group = q{};

my %remap_groups = ( '0.0' => 1,
                     '0.1' => 2,
                     '1.0' => 3,
                     '1.1' => 4,
                     '2.0' => 5,
                     '2.1' => 6, );

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (.+ \s) (\d\.\d) \s* \z/xms ) { 
        $front = $1;
        $group = $2;
        if ( exists $remap_groups{$group} ) { 
            $group = $remap_groups{$group};
        }
        print "$front$group\n";
    }
}

