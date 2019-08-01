#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;

# Given file names of input reads, I want output script lines sort of like this:
#    salmon --no-version-check quant --threads 6 --seqBias --gcBias --validateMappings --libType A \
#        --index TAIR10_cdna_20101214_updated_w_GFP.salmon_index.2018.01.09 \
#        --geneMap TAIR10_cdna_20101214_updated_w_GFP.tx2gene.tsv.txt \
#        --output [output name] \
#        --unmatedReads ../filtered_reads1/replicate_3nt.exact_1/Col_WT_rep1.trim_exact_3nt.fq.gz --numBootstraps 200 ;

my $header = '#!/bin/bash' 
            . "\n\n" 
            . ". \"/home/bioinformatics/anaconda2/etc/profile.d/conda.sh\" ;\n\n" 
            . "conda activate salmon_0.14.1 ;\n\n"
            ;

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
        $output = "$stem.salmon_2019.07.31.01";
    }
    else {
        die "Cannot parse basename $basename (from input $input)\n"
    }

    $output = safename($output);

    print $header if $header;
    $header = q{};

    print "salmon --no-version-check quant --threads 6 --seqBias --gcBias --validateMappings --libType A ",
          "--index /home/bioinformatics/sepal_rnaseq_dec2016/salmon_cDNA_w_UTRs/TAIR10_cdna_20101214_updated_w_GFP.salmon_index.2019.07.31 ",
          "--geneMap /home/bioinformatics/sepal_rnaseq_dec2016/salmon_cDNA_w_UTRs/TAIR10_cdna_20101214_updated_w_GFP.tx2gene.tsv.txt ",
          "--output $output --unmatedReads $input --numBootstraps 200 ;",
          "\n",
          ;
}

print "\nconda deactivate ;\n\n";

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

