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

        my $job      = 'job_kallisto_sepals_' . $read_tag . "_2016.01.25.$j";
        my $job_file = "$job.sh";
        $job_file    = safename($job_file);

        open my $JOB, '>', $job_file;

        print $JOB '#!/bin/bash -login', "\n";
        print $JOB '#PBS -l walltime=006:00:00', "\n";
        print $JOB '#PBS -l nodes=1:ppn=8', "\n";
        print $JOB '#PBS -l mem=32gb', "\n";
        print $JOB "#PBS -N $job.sh\n";
        print $JOB '#PBS -q main', "\n";
        print $JOB '#PBS -M ems394@cornell.edu', "\n";
        print $JOB '#PBS -m abe', "\n";
        print $JOB '#PBS -A ged', "\n";
        print $JOB '#PBS -r n', "\n";
        print $JOB '#PBS -V', "\n";
        print $JOB 'cd /mnt/ls15/scratch/users/emsch/work_rsync/2015/adrienne/rsem_2016 ;', "\n";
        print $JOB 'kallisto quant';
        print $JOB ' -i /mnt/home/emsch/work/2015/adrienne/rsem_2015/indices/Athal_sepal_seqs_2015.07.26.kallisto'; 
        print $JOB ' -t 8 --single';
        if ( $size > 0 ) {
            print $JOB " -l $size -s $sd";
        }
        if ( $size <= 0 ) {
            die "Without defined mean insert length, kallisto will fail in single-end mode for $reads\n";
        }
        print $JOB ' -b 100 -o kallisto_', $read_tag, "_2016.01.25.$j";
        print $JOB ' --pseudobam';
        print $JOB " $reads";
        print $JOB ' 1>kallisto_', $read_tag, "_2016.01.25.$j.bam";
        print $JOB ' 2>kallisto_', $read_tag, "_2016.01.25.$j.err ;\n";
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

