#!/usr/bin/env perl

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

while (my $input = <>) { 
    chomp $input;
    my @vals     = split /\t/, $input;
    my @new_vals = ();
    foreach my $val (@vals) { 
        if ( ( looks_like_number($val) ) and ( $val == 0 ) ) { 
            $val = 0.001;
        }
        push @new_vals, $val;
    }
    my $output = join "\t", @new_vals;
    print "$output\n";
}

