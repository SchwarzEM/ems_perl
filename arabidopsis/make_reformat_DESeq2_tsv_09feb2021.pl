#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = '#!/bin/bash';

while (my $input = <>) {
    chomp $input; 
    # sample input: 
    # 8hrMock.vs.0hr_DESeq2_all.data_2021.02.08.01.orig.tsv.txt
    if ( $input =~ /\A ((\S+) _DESeq2_all\.data (?:_|\.) \d+\.\d+\.\d+\.\d+) \.orig \.tsv \.txt \z/xms ) { 
        my $file_stem = $1;
        my $geno_stem = $2;

        # Print header exactly once, at the top of the output:
        print "$header\n\n" if $header;
        $header = q{};

        print "cat $input | ";
        print 'perl -ne \' s/\A["]["]/Gene/; s/["]//g; ';
        print 's/\blog2FoldChange\b/', $geno_stem, '.log2FoldChange/g; ';
        print 's/\bpadj\b/', $geno_stem, '.padj/g; ';
        print "print; ' | cut -f 1,3,7 > $file_stem.tsv.txt ;\n"; 
        print "\n";
    }
    else {
        die "Cannot parse input: $input\n";
    }
}
