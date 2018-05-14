#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $i = 0;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \s+ (\d+) \s+ (\d+) \z/xms ) {
        my $reads  = $1;
        my $size   = $2;
        my $sd     = $3;

        if (! -e $reads) {
            die "Cannot find readfile $reads\n";
        }

        my $read_tag = q{};
        if ( $reads =~ /\A \S+ \/ ([^\/\s]+) \. (?: fq|fastq) \z/xms ) {
            $read_tag = $1;
        }
        else {
            die "Can't parse readfile name $reads\n";
        }

        $i++;
        my $j = sprintf "%02i", $i;

        my $job      = 'job_rsem_sepals_' . $read_tag . "_2016.01.21.$j";
        my $job_file = "$job.sh";
        $job_file    = safename($job_file);

        open my $JOB, '>', $job_file;

        print $JOB '#!/bin/bash -login', "\n";
        print $JOB '#PBS -l walltime=006:00:00', "\n";
        print $JOB '#PBS -l nodes=1:ppn=8', "\n";
        print $JOB '#PBS -l mem=32gb', "\n";
        print $JOB "#PBS -N $job\n";
        print $JOB '#PBS -q main', "\n";
        print $JOB '#PBS -M ems394@cornell.edu', "\n";
        print $JOB '#PBS -m abe', "\n";
        print $JOB '#PBS -A ged', "\n";
        print $JOB '#PBS -r n', "\n";
        print $JOB '#PBS -V', "\n";
        print $JOB 'cd /mnt/ls15/scratch/users/emsch/work_rsync/2015/adrienne/rsem_2016 ;', "\n";
        print $JOB 'PATH="$PATH:/mnt/home/emsch/src/bowtie2-2.2.5/bin" ;', "\n";
        print $JOB 'PATH="$PATH:/mnt/home/emsch/src/rsem-1.2.21/bin" ;', "\n";
        print $JOB 'module load SAMTools/1.2 ;', "\n";
        print $JOB 'rsem-calculate-expression --bowtie2 -p 8 --calc-pme --phred33-quals --calc-ci --ci-credibility-level 0.99 --ci-memory 30000';
        if ( $size > 0 ) {
            print $JOB " --fragment-length-mean $size --fragment-length-sd $sd";
        }
        print $JOB " $reads";
        print $JOB ' /mnt/ls15/scratch/users/emsch/work_rsync/2015/adrienne/rsem_2016/indices/Athal_sepal_seqs_2015.07.26.rsem';
        print $JOB ' rsem_', $read_tag, "_2016.01.21.$j ;\n";

        close $JOB;

    }
    else {
        die "Can't parse input: $input\n";
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

