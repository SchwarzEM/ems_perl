#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    while ( $input =~ /\A \S+? ([N]+) (.*) \z/xmsg ) { 
        my $n_block   = $1;
        my $remainder = $2;

        my $n_block_count = ($n_block =~ s/N/N/g);
        print "$n_block_count\n";

        $input = $remainder;
    }
}

