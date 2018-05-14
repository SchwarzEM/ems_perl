#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A [#] \s+ (\S+) \s+ vs\. \s+ (\S+) \s* \z/xms ) { 
        my $top    = $1;
        my $bottom = $2;
        print "$input\n";
        print q{> comp <- exactTest(edgeR_geno_comp, pair=c("}, $bottom, q{","}, $top, q{"))}, "\n";
        print "> deg <- topTags(comp, n=Inf, p.value=0.1)\n";
        print "> summary(decideTestsDGE(comp, p.value=0.1))\n";
        print q{> write.csv(as.data.frame(deg), file="}, $top, q{_}, $bottom, q{_edgeR_exactTest_padj0.1_2016.05.20.csv")}, "\n";
        print "\n";
    }
    else { 
        die "Cannot parse input: $input\n";
    }
}


