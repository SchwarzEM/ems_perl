#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = "Gene\tSymbols\tAnnotation\n";

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        my $gene    = q{};
        my $symbols = q{};
        my $annot   = q{};
        if ( $input =~ /\A > (AT\S+G\d+) \. \d+ \s+ \| \s+ Symbols: \s (.*) \s+ \| \s+ (.+?) \s+ \|  /xms ) { 
            $gene    = $1;
            $symbols = $2;
            $annot   = $3;
            print $header if $header;
            $header = q{};
            print "$gene\t$symbols\t$annot\n";
        }
        else {
            die "Can't parse header: $input\n";
        }
    }
}
        
