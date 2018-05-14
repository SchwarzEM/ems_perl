#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A \S+ \/ (jj\d+) \.bowtie2\.WBcel235\.bt2\.sorted\.bam \z/xms ) {
        my $allele   = $1;

        my $outfile = "$allele.bowtie2.WBcel235.bt2.sorted.dedup.bam";
        $outfile    = safename($outfile);

        my $job1 = 'job_samtools_' . $allele . '_dedup_2016.10.03.01.sh';

        $job1 = safename($job1);

        open my $JOB1, '>', $job1 ;

        print $JOB1 '#!/bin/bash -login', "\n";
        print $JOB1 '#PBS -l walltime=003:59:00', "\n";
        print $JOB1 '#PBS -l nodes=1:ppn=8', "\n";
        print $JOB1 '#PBS -l mem=32gb', "\n";
        print $JOB1 "#PBS -N $job1\n";
        print $JOB1 '#PBS -q main', "\n";
        print $JOB1 '#PBS -M ems394@cornell.edu', "\n";
        print $JOB1 '#PBS -m abe', "\n";
        print $JOB1 '#PBS -A ged', "\n";
        print $JOB1 '#PBS -r n', "\n";
        print $JOB1 '#PBS -V', "\n";
        print $JOB1 "cd /mnt/ls15/scratch/users/emsch/kelly_mut_seqs ;\n";
        print $JOB1 "module load SAMTools/1.3.1 ;\n";
        print $JOB1 "samtools rmdup -s $input $outfile ;\n";

        close $JOB1;
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

