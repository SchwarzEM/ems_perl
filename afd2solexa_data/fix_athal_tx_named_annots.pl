#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = "Gene\tAnnotation";

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (AT (?: \d|M|C) G \d+) \. \d \t (.*) \z/xms ) { 
        my $gene  = $1;
        my $annot = $2;
        print "$header\n" if $header;
        $header = q{};
        print "$gene\t$annot\n";
    }
    else { 
        if ( $input !~ /\A Gene \t/xms ) { 
            warn "Unparsed:  $input\n";
        }
    }
}

