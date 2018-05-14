#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;

# Given file names of input reads, I want output script lines sort of like this:
#    salmon --no-version-check quant -p 6 --seqBias -l A -i TAIR10_cds_20101214_updated_w_GFP.salmon_index.2016.12.12 \
#        -g TAIR10_cds_20101214_updated_w_GFP.tx2gene.tsv.txt -o [output name] -r ../filtered_reads1/replicate_3nt.exact_1/Col_WT_rep1.trim_exact_3nt.fq.gz --numBootstraps 100 ;

my $header = '#!/bin/bash' . "\n\n";

while (my $input = <>) {
    chomp $input;

    my $output = q{};

    # Sample input names:
    # /home/bioinformatics/sepal_rnaseq_dec2016/filtered_reads1/*/*.fq.gz
    # [or]
    # /home/bioinformatics/sepal_rnaseq_dec2016/filtered_reads2b/*.gz

    if (! -r $input) {
        die "Cannot read input file: $input\n";
    }

    my $basename = basename($input);
    if ( $basename =~ /\A (\S+) \.fq\.gz \z/xms ) { 
        my $stem = $1;
        $output = "$stem.salmon_08jul2017";
    }
    else {
        die "Cannot parse basename $basename (from input $input)\n"
    }

    $output = safename($output);

    print $header if $header;
    $header = q{};

    print "salmon --no-version-check quant -p 6 --seqBias -l A",
          " -i /home/bioinformatics/sepal_rnaseq_dec2016/salmon/TAIR10_cds_20101214_updated_w_GFP.salmon_index.2016.12.12",
          " -g /home/bioinformatics/sepal_rnaseq_dec2016/salmon/TAIR10_cds_20101214_updated_w_GFP.tx2gene.tsv.txt",
          " -o $output -r $input --numBootstraps 100 --writeUnmappedNames ;",
          "\n",
          ;
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

