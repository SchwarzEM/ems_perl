#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Getopt::Long;
use File::Basename;

# Given file names of input reads, I want output script lines like this:
# 
#     salmon --no-version-check quant --threads 4 --validateMappings --seqBias --gcBias --libType A \
#     --index /home/bioinformatics/sepal_rnaseq_dec2016/salmon_cDNA_w_UTRs/TAIR10_cdna_20101214_updated_w_GFP_gentrome_index \
#     --geneMap /home/bioinformatics/sepal_rnaseq_dec2016/salmon_cDNA_w_UTRs/TAIR10_cdna_20101214_updated_w_GFP.tx2gene.tsv.txt \
#     --unmatedReads 8_hr_11-24-20_fastp.filt/12048_11190_126769_H3G7KBGXH_1uM3-8hrs_GATCAG_R1.fastp_filt1.fq.gz \
#     --output 8_hr_11-24-20_salmon/12048_11190_126769_H3G7KBGXH_1uM3-8hrs_GATCAG_R1.fastp_filt1.salmon ;

my @infiles = ();
my $outdir  = q{};
my $threads = q{};
my $index   = q{};
my $geneMap = q{};

my $header          = '#!/bin/bash' . "\n\n";
my $footer          = q{};

my $help;

GetOptions ( 'infiles=s{,}' => \@infiles,
             'outdir=s'     => \$outdir,
             'threads=i'    => \$threads,
             'index=s'      => \$index,
             'geneMap=s'    => \$geneMap,
             'help'         => \$help,   );

if ( $help or (! @infiles ) or (! $threads ) or (! $index ) or (! $geneMap ) ) { 
    die "Format: make_salmon_job_12-28-20.pl\n",
        "    --infile|-i    <input stream/files>\n",
        "    --outdir|-o    <directory for all output files>\n",
        "    --threads|-t   <threads argument for salmon>\n",
        "    --index|-i     <index argument for salmon>\n",
        "    --geneMap|-g   <geneMap argument for salmon>\n",
        "    --help|-h      [print this message]\n",
        ;
}

if (! -d $index ) {
    die "Do not observe a salmon index directory called: $index\n";
}

if (! -r $geneMap ) {
    die "Cannot read geneMap file at: $geneMap\n";
}

foreach my $input (@infiles) { 
    my $output = q{};
    my $stem   = q{};
    my $suffix = q{};

    # Sample input name:
    # 8_hr_11-24-20_fastp.filt/12048_11190_126769_H3G7KBGXH_1uM3-8hrs_GATCAG_R1.fastp_filt1.fq.gz

    if (! -r $input) {
        die "Cannot read input file: $input\n";
    }

    my $basename = basename($input);
    if ( $basename =~ /\A (\S+) \.fq (?:\.gz){0,1} \z/xms ) { 
        $stem   = $1;
        $output = "$stem.salmon";
    }
    else {
        die "Cannot parse basename $basename (from input $input)\n"
    }

    print $header if $header;
    $header = q{};
    $footer = "\n";

    print "salmon --no-version-check quant";
    print " --threads $threads --validateMappings --seqBias --gcBias --libType A";
    print " --index $index";
    print " --geneMap $geneMap";
    print " --unmatedReads $input";
    print " --output $outdir/$output";
    print " ;\n"; 
}

print $footer if $footer;
