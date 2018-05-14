#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Getopt::Long;
use File::Basename;
use Scalar::Util qw(looks_like_number);

my $salmon_index   = q{};
my $salmon_tx2gene = q{};
my @read_files     = ();
my $bootstraps     = q{};
my $output_tag     = q{};

my $help;

GetOptions ( 
    'index=s'         => \$salmon_index,
    'tx2gene=s'       => \$salmon_tx2gene,
    'read_files=s{,}' => \@read_files,
    'bootstraps=i'    => \$bootstraps,
    'output_tag=s'    => \$output_tag,
    'help'            => \$help, 
);

if ( $help 
     or (! -r $salmon_index ) 
     or (! -r $salmon_tx2gene ) 
     or (! @read_files ) 
     or (! looks_like_number($bootstraps) ) 
     or ( $bootstraps != int($bootstraps) )
     or ( $bootstraps < 1 )
   ) {
    die "Format: make_salmon_job_12jan2018.pl\n",
        "    --index|-i       [Salmon transcript index]\n",
        "    --tx2gene|-t     [Salmon transcript-to-gene TSV table]\n",
        "    --reads|-r       [input read files]\n",
        "    --bootstraps|-b  [number of bootstraps (positive integer)]\n",
        "    --output_tag|-o  [suffix to mark output, e.g, \"cds.12jan2018\"]\n",
        "    --help|-h        [print this message]\n",
        ;
}

my $header = '#!/bin/bash' . "\n\n";

foreach my $read_file (@read_files) {
    chomp $read_file;
    if (! -r $read_file ) {
        die "Cannot read input RNA-seq file: $read_file\n";
    }

    my $output = q{};

    # Sample input names:
    # /home/bioinformatics/sepal_rnaseq_dec2016/filtered_reads1/*/*.fq.gz
    # [or]
    # /home/bioinformatics/sepal_rnaseq_dec2016/filtered_reads2b/*.gz

    my $basename = basename($read_file);
    if ( $basename =~ /\A (\S+) \.fq\.gz \z/xms ) { 
        my $stem = $1;
        $output = "$stem.salmon.$output_tag";
    }
    else {
        die "Cannot parse basename $basename (from input read file $read_file)\n"
    }

    $output = safename($output);

    # Do this exactly once at the start of the outputs:
    print $header if $header;
    $header = q{};

    print "salmon",
          " --no-version-check",
          " quant",
          " -p 6",
          " -i $salmon_index",
          " -g $salmon_tx2gene",
          " -l A",
          " -r $read_file",
          " -o $output",
          " --seqBias --gcBias",
          " --numBootstraps $bootstraps",
          " ;",
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

