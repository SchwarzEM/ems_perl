#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = '#!/bin/bash';

while (my $input = <>) {
    chomp $input; 
    # sample input: 
    # 12dpi_noDEX.vs.12.D_edgeR_exactTest_all.data_2019.06.01.01.csv
    if ( $input =~ /\A ((\S+) _edgeR_exactTest_all\.data (?:_|\.) \d+\.\d+\.\d+\.\d+)\.csv \z/xms ) { 
        my $file_stem = $1;
        my $geno_stem = $2;

        # Print header exactly once, at the top of the output:
        print "$header\n\n" if $header;
        $header = q{};

        print "cat $input | ";
        print 'perl -ne \' s/\A["]["]/Gene/; s/[,]/\t/g; s/["]//g; ';
        print 's/\blogFC\b/', $geno_stem, '.logFC/g; ';
        print 's/\blogCPM\b/', $geno_stem, '.logCPM/g; ';
        print 's/\bPValue\b/', $geno_stem, '.PValue/g; ';
        print 's/\bFDR\b/', $geno_stem, '.FDR/g; ';
        print "print; ' | cut -f 1-2,5 > $file_stem.tsv.txt ;\n"; 
        print "\n";
    }
    else {
        die "Cannot parse input: $input\n";
    }
}
