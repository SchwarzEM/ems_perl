#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;

    # Sample input line:
    # VCFs/CB4856.bowtie2.WBcel235.bt2.sorted.dedup.vcf

    if ( $input =~ /\A \S+ \/ ( ([A-Za-z0-9]+) \. \S+) \.vcf \z/xms ) { 
        my $filestem = $1;
        my $allele   = $2;
        my $out_pref = "$allele.bowtie2.WBcel235";

        my $job1 = 'job_bcftools.convert_' . $allele . '_2016.10.05.01.sh';
        my $job2 = 'job_bcftools.index_' . $allele . '_2016.10.05.01.sh';

        $job1 = safename($job1);
        $job2 = safename($job2);

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
        print $JOB1 "bcftools convert --threads 8 --output-type b --output $filestem.bcf $input ;\n";
        print $JOB1 "qsub $job2 ;\n";

        close $JOB1;
        open my $JOB2, '>', $job2;

        print $JOB2 '#!/bin/bash -login', "\n";
        print $JOB2 '#PBS -l walltime=003:59:00', "\n";
        print $JOB2 '#PBS -l nodes=1:ppn=1', "\n";
        print $JOB2 '#PBS -l mem=32gb', "\n";
        print $JOB2 "#PBS -N $job2\n";
        print $JOB2 '#PBS -q main', "\n";
        print $JOB2 '#PBS -M ems394@cornell.edu', "\n";
        print $JOB2 '#PBS -m abe', "\n";
        print $JOB2 '#PBS -A ged', "\n"; 
        print $JOB2 '#PBS -r n', "\n";
        print $JOB2 '#PBS -V', "\n";
        print $JOB2 "cd /mnt/ls15/scratch/users/emsch/kelly_mut_seqs ;\n";
        print $JOB2 "bcftools index $filestem.bcf ;\n";

        close $JOB2;
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

