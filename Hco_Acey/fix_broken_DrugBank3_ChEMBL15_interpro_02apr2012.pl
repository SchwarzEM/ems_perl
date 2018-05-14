#!/usr/bin/env perl

use strict;
use warnings;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A ([^\t]+ \t [^\t]+ \t [^\t]+ \t) (\d+) \z/xms ) { 
        my $text1 = $1;
        my $text2 = $2;
        print $text1, 'drugbank_target|', "$text2\n", ;
    }
    else { 
        print "$input\n";
    }
}

