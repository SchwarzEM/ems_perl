#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <> ) {
    chomp $input;
    # Sample input:
    # Scarp_III_RagTag  Liftoff mRNA    9096    12050   .       +       .       ID=transcript:L596_016927;Parent=gene:L596_016927;
    if ( $input =~ /\A \S+ \t \S+ \t mRNA \t \d+ \t \d+ \t [^\t]* \t [+-] \t [^\t]* \t ID=transcript: ([^;\s]+) [;] Parent=gene: ([^;\s]+) [;] /xms ) { 
         my $tx   = $1;
         my $gene = $2;
         print "$tx\t$gene\n";
    }
    elsif ( $input =~ /\A [^\t]* \t [^\t]* \t mRNA \t /xms ) {
        die "Cannot parse input mRNA line: $input\n";
    }
}
