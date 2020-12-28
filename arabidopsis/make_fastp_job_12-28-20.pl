#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Getopt::Long;
use File::Basename;

# Given file names of input reads, I want output script lines like this:
# 
#     fastp --thread 4 --dont_overwrite --length_required 75 --max_len1 84 \
#     --json 12048_11190_126761_H3G7KBGXH_Mock1-8hrs_ATCACG_R1.json \
#     --html 12048_11190_126761_H3G7KBGXH_Mock1-8hrs_ATCACG_R1.html \
#     --in1 ./filtered_reads/12048_11190_126761_H3G7KBGXH_Mock1-8hrs_ATCACG_R1.filt1.fq \
#     --out1 ./fastp/12048_11190_126761_H3G7KBGXH_Mock1-8hrs_ATCACG_R1.filt1.fastp.fq ;

my @infiles         = ();
my $outdir          = q{};
my $threads         = q{};
my $length_required = q{};
my $max_len1        = q{};

my $header          = '#!/bin/bash' . "\n\n";
my $footer          = q{};

my $help;

GetOptions ( 'infiles=s{,}'      => \@infiles,
             'outdir=s'          => \$outdir,
             'thread=i'          => \$threads,
             'length_required=i' => \$length_required,
             'max_len1=i'        => \$max_len1,
             'help'              => \$help,   );

if ( $help or (! @infiles ) or (! $threads ) or (! $length_required ) or (! $max_len1 ) ) { 
    die "Format: make_fastp_job_12-28-20.pl\n",
        "    --infile|-i            <input stream/files>\n",
        "    --outdir|-o            <directory for all output files>\n",
        "    --thread|-t            <thread argument, i.e., number of threads for fastp>\n",
        "    --length_required|-l   <length_required argument for fastp>\n",
        "    --max_len1|-m          <max_len1 argument for fastp>\n",
        "    --help|-h              [print this message]\n",
        ;
}

foreach my $input (@infiles) { 
    my $output = q{};
    my $stem   = q{};
    my $suffix = q{};

    # Sample input name:
    # /home/bioinformatics/Vijaya/ATML1_induction_RNA-seq/filtered_reads/12048_11190_126761_H3G7KBGXH_Mock1-8hrs_ATCACG_R1.filt1.fq

    if (! -r $input) {
        die "Cannot read input file: $input\n";
    }

    my $basename = basename($input);
    if ( $basename =~ /\A (\S+) \.filt1 \. (fq (?:\.gz){0,1} ) \z/xms ) { 
        $stem   = $1;
        $suffix = $2;
        $output = "$stem.fastp_filt1.$suffix";
    }
    else {
        die "Cannot parse basename $basename (from input $input)\n"
    }

    print $header if $header;
    $header = q{};
    $footer = "\n";

    print "fastp --thread $threads --dont_overwrite";
    print " --length_required $length_required --max_len1 $max_len1";
    print " --json $outdir/$stem.json";
    print " --html $outdir/$stem.html";
    print " --in1 $input";
    print " --out1 $outdir/$output";
    print " ;\n"; 
}

print $footer if $footer;
