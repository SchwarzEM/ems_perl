#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = '#!/bin/bash';

while (my $input = <>) {
    chomp $input; 
    # sample input: 
    # Hco_fem_worm.vs.Hco_fem_gut_edgeR_exactTest_all.data_2018.09.19.01.csv
    if ( $input =~ /\A ((\S+) _edgeR_exactTest_all\.data_\S+)\.csv \z/xms ) { 
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
        print "print; ' > $file_stem.tsv.txt ;\n"; 
        print "\n";
    }
    else {
        die "Cannot parse input: $input\n";
    }
}
