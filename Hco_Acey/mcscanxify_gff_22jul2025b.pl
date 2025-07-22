#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while ( my $input = <> ) {
    chomp $input;
    # Sample input:
    # Sherm_X   AUGUSTUS        gene    15650   20230   .       -       .       ID=Sherm_X.g16154;
    if ( $input =~ /\A (\S+) \t [^\t]* \t gene \t (\d+) \t (\d+) \t .* \t ID[=] ([^;\s]+) [;] /xms ) { 
        my $chr      = $1;
        my $start_nt = $2;
        my $end_nt   = $3;
        my $gene     = $4;
        print "$chr\t$gene\t$start_nt\t$end_nt\n";
    }
}
