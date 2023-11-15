#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    # Sample input:
    # Sherm_2022.07.24.06     9049    9549    .       .       +       AUGUSTUS        gene    .       ID=Sherm_2022.07.24.06.g2;
    chomp $input;
    if ( $input =~ /\A ([^\t]* \t [^\t]* \t [^\t]*) \t \. \t (.* \t ID = (\S+) [;])\z/xms ) {
        my $front = $1;
        my $back  = $2;
        my $gene  = $3;
        print "$front\t$gene\t$back\n";
    }
    else {
        die "Cannot parse: $input\n";
    }
}
