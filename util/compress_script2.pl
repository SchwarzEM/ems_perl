#!/usr/bin/env perl

# compress_script2.pl -- Erich Schwarz <emsch@caltech.edu>, 10/6/2011.
# Purpose: try getting more reliable compression of shell scripts.

use strict;
use warnings;

my $output_text = q{};

while (my $input = <>) { 
    if ( $input =~ / \A (.* \s) \\ \s* \n /xms ) { 
        $input = $1;
    }
    $output_text .= $input;
}

$output_text =~ s/[ ]+/ /g;

print "$output_text";


