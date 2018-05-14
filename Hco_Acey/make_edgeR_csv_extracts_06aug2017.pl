#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header  = '#!/bin/bash' . "\n\n";
my $outfile = q{};

while (my $infile = <>) {
    chomp $infile;
    if ( $infile =~ /\A \S+ \/ ((\S+?) _edgeR\S+) \.csv \z/xms ) { 
        my $stem = $1;
        my $tag  = $2;

        my $tag_logFC = $tag . '_log2FC';
        my $tag_FDR   = $tag . '_FDR';

        $outfile = "$stem.tsv.txt";
        $outfile    = safename($outfile);

        print $header if $header;
        $header = q{};
        
        print "cat $infile | ";
        print "perl -ne ";
        print q{' s/[,]/\t/g; s/["]{2}/Gene/g; s/["]//g; s/logFC/};
        print "$tag_logFC";
        print q{/g; s/FDR/};
        print "$tag_FDR";
        print q{/g; print; '};
        print " | cut -f 1-2,5 > $outfile ;\n";
    }
    else {
        die "Cannot parse input file name: $infile\n";
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

