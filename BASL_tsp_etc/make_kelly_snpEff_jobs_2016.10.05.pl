#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;

    # sample input that I want:
    # select_VCFs/jj103.not.CB4856.bowtie2.WBcel235.bt2.sorted.dedup.vcf

    if ( $input =~ /\A \S+ \/ (\S+)\.not\.CB4856\.bowtie2\.WBcel235\.bt2\.sorted\.dedup\.vcf \z/xms ) {
        my $genotype = $1;

        my $outfile1 = "$genotype.not.CB4856.bowtie2.WBcel235.bt2.sorted.dedup.ann_eff.vcf";
        $outfile1    = safename($outfile1);

        my $err1 = "$genotype.not.CB4856.bowtie2.WBcel235.bt2.sorted.dedup.ann_eff.err";
        $err1 = safename($err1);

        my $outfile2 = "$genotype.not.CB4856.bowtie2.WBcel235.bt2.sorted.dedup.ann_classic.vcf";
        $outfile2    = safename($outfile2);
        
        my $err2 = "$genotype.not.CB4856.bowtie2.WBcel235.bt2.sorted.dedup.ann_classic.err";
        $err2 = safename($err2);

        my $csvstats = "$genotype.not.CB4856.bowtie2.WBcel235.bt2.sorted.dedup.ann_eff.stats.csv";
        $csvstats    = safename($csvstats);

        my $stats_summ = "$genotype.not.CB4856.bowtie2.WBcel235.bt2.sorted.dedup.ann_eff.summary.html";
        $stats_summ    = safename($stats_summ);

        my $jobfile = 'job_snpEff_' . $genotype . '.not.CB4856.bowtie2.WBcel235_2016.10.04.01.sh';
        $jobfile    = safename($jobfile);

        open my $JOB, '>', $jobfile;

        print $JOB '#!/bin/bash -login', "\n";
        print $JOB '#PBS -l walltime=003:59:00', "\n";
        print $JOB '#PBS -l nodes=1:ppn=1', "\n";
        print $JOB '#PBS -l mem=32gb', "\n";
        print $JOB "#PBS -N $jobfile\n";
        print $JOB '#PBS -q main', "\n";
        print $JOB '#PBS -M ems394@cornell.edu', "\n";
        print $JOB '#PBS -m abe', "\n";
        print $JOB '#PBS -A ged', "\n";
        print $JOB '#PBS -r n', "\n";
        print $JOB '#PBS -V', "\n";
        print $JOB 'cd /mnt/ls15/scratch/users/emsch/kelly_mut_seqs ;', "\n";
        print $JOB 'module load Java/jdk1.8.0 ;', "\n";
        print $JOB "java -Xmx31g -jar /mnt/home/emsch/src/snpEff_04jul2016/snpEff.jar ";
        print $JOB "eff -lof -csvStats $csvstats -htmlStats $stats_summ WBcel235.82 $input 1>", "$outfile1 2>", "$err1 ;\n";
        print $JOB "java -Xmx31g -jar /mnt/home/emsch/src/snpEff_04jul2016/snpEff.jar ";
        print $JOB "eff -lof -classic WBcel235.82 $input 1>", "$outfile2 2>", "$err2 ;\n";

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

