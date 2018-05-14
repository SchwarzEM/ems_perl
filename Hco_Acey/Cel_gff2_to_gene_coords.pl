#!/usr/bin/env perl

use strict;
use warnings;

# CHROMOSOME_II	gene	gene	9197458	9207148	.	+	.	Gene "WBGene00002299" ; Position "1.09455" ; Locus "let-23"
# CHROMOSOME_III	gene	gene	7528603	7536567	.	-	.	Gene "WBGene00003024" ; Position "-0.674053" ; Locus "lin-39"

my $header = "Cel_Gene\tCel_coords\n";

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A CHROMOSOME_(\S+) \t gene \t gene \t (\d+) \t (\d+) \t [^\t]* \t (\S+) \t [^\t]* \t Gene [ ] \" (WBGene\d+) \" /xms ) {
        my $chr        = $1;
        my $five_p_nt  = $2;
        my $three_p_nt = $3;
        my $ori        = $4;
        my $wbgene     = $5;
        print $header if $header;
        $header = q{};
        print $wbgene, "\t", $chr, q{:}, $five_p_nt, q{-}, $three_p_nt, " [$ori]\n", ;
    }
}

