#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A \S+ \/ (\S+)\.bowtie2\.WBcel235\.bt2\.sorted\.dedup\.bam \z/xms ) { 
        my $genotype = $1;

        my $outfile  = "$genotype.bowtie2.WBcel235.bt2.sorted.dedup.vcf";
        $outfile     = safename($outfile);

        my $jobfile = 'job_samtools_' . $genotype . '_vcf_2016.10.04.01.sh';
        $jobfile    = safename($jobfile);

        open my $JOB, '>', $jobfile;

        print $JOB '#!/bin/bash -login', "\n";
        print $JOB '#PBS -l walltime=003:59:00', "\n";
        print $JOB '#PBS -l nodes=1:ppn=1', "\n";
        print $JOB '#PBS -l mem=32gb', "\n";
        print $JOB "#PBS -N $jobfile \n";
        print $JOB '#PBS -q main', "\n";
        print $JOB '#PBS -M ems394@cornell.edu', "\n";
        print $JOB '#PBS -m abe', "\n";
        print $JOB '#PBS -A ged', "\n";
        print $JOB '#PBS -r n', "\n";
        print $JOB '#PBS -V', "\n";
        print $JOB 'cd /mnt/ls15/scratch/users/emsch/kelly_mut_seqs ;', "\n";
        print $JOB 'module load SAMTools/1.2 ;', "\n";
        print $JOB 'samtools mpileup -vu -f Caenorhabditis_elegans.WBcel235.dna.toplevel.fa';
        print $JOB " -o $outfile $input ;\n";

        close $JOB;
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


