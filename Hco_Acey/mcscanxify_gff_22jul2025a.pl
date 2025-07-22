#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while ( my $input = <> ) {
    chomp $input;
    # Sample input:
    # Scarp_III_RagTag        Liftoff gene    9096    12050   .       +       .       ID=gene:L596_016927; [...]
    if ( $input =~ /\A (\S+) \t [^\t]* \t gene \t (\d+) \t (\d+) \t .* \t ID[=]gene: ([^;\s]+) [;] /xms ) { 
        my $chr      = $1;
        my $start_nt = $2;
        my $end_nt   = $3;
        my $gene     = $4;
        print "$chr\t$gene\t$start_nt\t$end_nt\n";
    }
}
