#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = '#!/bin/bash';

while (my $input = <>) {
    chomp $input; 
    # sample input: 
    # WT.vs.vos2_DESeq2_2016.05.16.csv
    if ( $input =~ /\A ((\S+) _DESeq2 (?:_|\.) \d+\.\d+\.\d+ (?:\.\d+){0,1} )\.csv \z/xms ) { 
        my $file_stem = $1;
        my $geno_stem = $2;

        # Print header exactly once, at the top of the output:
        print "$header\n\n" if $header;
        $header = q{};

        print "cat $input | ";
        print 'perl -ne \' s/\A["]["]/Gene/; s/[,]/\t/g; s/["]//g; ';
        print 's/\bbaseMean\b/', $geno_stem, '.baseMean/g; ';
        print 's/\blog2FoldChange\b/', $geno_stem, '.log2FoldChange/g; ';
        print 's/\blfcSE\b/', $geno_stem, '.lfcSE/g; ';
        print 's/\bstat\b/', $geno_stem, '.stat/g; ';
        print 's/\bpvalue\b/', $geno_stem, '.pvalue/g; ';
        print 's/\bpadj\b/', $geno_stem, '.padj/g; ';
        print "print; ' > $file_stem.tsv.txt ;\n"; 
        print "\n";
    }
    else {
        die "Cannot parse input: $input\n";
    }
}

