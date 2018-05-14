#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = '#!/bin/bash' . "\n\n";

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A \S+ \/ (\S+) .not.CB4856.bowtie2.WBcel235.bt2.sorted.dedup.ann_classic.vcf \z/xms ) { 
        my $genotype = $1;

        my $summary = "$genotype.key_gene_summary_06oct2016.txt";
        $summary    = safename($summary);

        print $header if $header;
        $header = q{};

        print "cat $input | grep -f key_genes.txt | grep -f key_flags.txt > $summary ;\n";
    }
    else { 
        die "Cannot parse input file $input\n";
    }
}

print "\n" if (! $header);

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


