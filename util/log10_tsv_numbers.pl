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
            if ( $input_val <= 0 ) { 
                die "Can't take log10 of a non-positive, non-zero number in: $input.\n";
            }
            $input_val = ( log($input_val) / log(10) );
        }
        push @output_vals, $input_val;
    }
    my $output = join "\t", @input_vals;
    print "$output\n";
}

