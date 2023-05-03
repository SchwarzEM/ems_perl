#!/usr/bin/env perl

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

while (my $input = <>) { 
    chomp $input;
    my @input_vals  = split "\t", $input;
    my @output_vals = ();
    foreach my $input_val (@input_vals) { 
        if ( looks_like_number($input_val) ) { 
            $input_val = ($input_val + 0.5);
            $input_val = int $input_val;
        }
        push @output_vals, $input_val;
    }
    my $output = join "\t", @input_vals;
    print "$output\n";
}

