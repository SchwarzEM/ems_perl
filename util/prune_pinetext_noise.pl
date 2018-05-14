#!/usr/bin/env perl

# prune_pinetext_noise.pl -- Erich Schwarz <emsch@its.caltech.edu>, 7/6/2010.
# Purpose: strip out junk from ends of lines in PINE-mail plaintext files.

use strict;
use warnings;

my $text = q{};

while (my $input = <>) { 
    chomp $input;
    $text = q{};
    if ( $input =~ /\A (.*) =20[=]* \s* \z /xms ) { 
        $text = $1;
        $input = $text;
    }
    print "$input\n";
}
    
