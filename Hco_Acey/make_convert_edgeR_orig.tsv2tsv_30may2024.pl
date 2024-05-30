#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = '#!/bin/bash';

while (my $input = <>) {
    chomp $input; 
    # sample inputs: 
    # 01/VT_7d_r_Fem.vs.VT_7d_m_Fem_edgeR_exactTest_2024.05.27.01.orig.tsv.txt
    if ( $input =~ /\A ((\S+) _edgeR_exactTest_ \d+\.\d+\.\d+\.\d+)\.orig\.tsv\.txt \z/xms ) { 
        my $file_stem = $1;
        my $geno_stem = $2;

        my $outfile = "$file_stem.tsv.txt";
        $outfile    = safename($outfile);

        # Print header exactly once, at the top of the output:
        print "$header\n\n" if $header;
        $header = q{};

        print "cat $input | ";
        print 'perl -ne \' s/\A["]["]/Gene/; s/["]//g; ';   # omit 's/[,]/\t/g; '
        print 's/\blogFC\b/', $geno_stem, '.logFC/g; ';
        print 's/\blogCPM\b/', $geno_stem, '.logCPM/g; ';
        print 's/\bPValue\b/', $geno_stem, '.PValue/g; ';
        print 's/\bFDR\b/', $geno_stem, '.FDR/g; ';
        print "print; ' > $outfile ;\n"; 
        print "\n";
    }
    else {
        die "Cannot parse input: $input\n";
    }
}

sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}


