#!/usr/bin/env perl

# compress_script.pl -- Erich Schwarz <emsch@caltech.edu>, 10/6/2011.
# Purpose: compress shell scripts from human-readable multi-line form (which bash doesn't seem to parse well) to less-readaable single-line form (usable for bash).

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


