#!/usr/bin/perl

use strict;
use warnings;
my %genes = ();

while (my $input = <>) {
    if ( $input !~ / :wp /xms ) {
        if ( $input =~ /\A ([^\t]+ \s [^\t]+) /xms ) {  
            $genes{$1} = 1;
        }
    }
}

foreach my $gene (sort keys %genes) {
    print "$gene\n";
}

