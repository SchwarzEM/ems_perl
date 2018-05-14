#!/usr/bin/env perl

use strict;
use warnings;

# Sample input line, from an AUGUSTUS 2.6.1 run:
# Acey_2012.08.05_0001	AUGUSTUS	gene	1	7875	0.3	-	.	Acey_2012.08.05_0001.g1

my $header = "Gene\tCoordinates\n";

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \t AUGUSTUS \t gene \t (\d+) \t (\d+) \t [^\t]* \t (\S+) \t [^\t]* \t (\S+\.g\d+) \s* \z/xms ) { 
        my $sequence = $1;
        my $five_p_nt    = $2;
        my $three_p_nt    = $3;
        my $ori      = $4;
        my $gene     = $5;
        print $header if $header;
        $header = q{};
        print $gene, "\t", $sequence, q{:}, $five_p_nt, q{-}, $three_p_nt, " [$ori]\n", ;
    }
}


