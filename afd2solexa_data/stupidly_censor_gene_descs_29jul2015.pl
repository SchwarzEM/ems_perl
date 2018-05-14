#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %seen = ();

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t .* \z/xms ) { 
        my $gene = $1;
        if (! $seen{$gene} ) { 
            print "$input\n";
            $seen{$gene} = 1;
        }
    }
    else { 
        die "Can't stupidly parse description line; $input\n";
    }
}


