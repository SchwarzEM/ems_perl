#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;

    # sample input that I want:
    # VCFs/jj103.bowtie2.WBcel235.bt2.sorted.dedup.vcf

    # file that I want to compare the sample to:
    # VCFs/CB4856.bowtie2.WBcel235.bt2.sorted.dedup.vcf

    if ( $input =~ /\A \S+ \/ (\S+)\.bowtie2\.WBcel235\.bt2\.sorted\.dedup\.bcf \z/xms ) { 
        my $genotype = $1;

        if ( $genotype eq 'CB4856' ) {
            die "No point in comparing CB4856 to itself; input: $input\n";
        }

        my $outfile  = "$genotype.not.CB4856.bowtie2.WBcel235.bt2.sorted.dedup.vcf";
        $outfile     = safename($outfile);

        my $jobfile = 'job_bcftools.isec_' . $genotype . '.not.CB4856_2016.10.05.01.sh';
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
        print $JOB 'bcftools isec --complement --collapse indels --collapse snps --output-type v';
        print $JOB " --output $outfile $input BCFs/CB4856.bowtie2.WBcel235.bt2.sorted.dedup.bcf ;\n";

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

