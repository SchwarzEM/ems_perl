#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input_genome = <>) {
    chomp $input_genome;

    my $stem = q{};

    if ( $input_genome !~ /\A \S+\.WS244\.genomic\.fa \z/xms ) {
        die "Can't parse input genome file name: $input_genome\n";
    }

    if ( $input_genome =~ /\A (\S+\.WS244)\.genomic\.fa \z/xms ) {
        $stem = $1;
    }

    my $job1 = 'job_RScout_' . $stem . '_20jul2014a.sh';
    my $job2 = 'job_RScout_' . $stem . '_20jul2014b.sh';
    my $job3 = 'job_RScout_' . $stem . '_20jul2014c.sh';

    $job1 = safename($job1);
    $job2 = safename($job2);
    $job3 = safename($job3);

    open my $JOB1, '>', $job1;

    print $JOB1 '#!/bin/bash -login', "\n";
    print $JOB1 '#PBS -l walltime=001:00:00', "\n";
    print $JOB1 '#PBS -l nodes=1:ppn=1', "\n";
    print $JOB1 '#PBS -l mem=25gb', "\n";
    print $JOB1 '#PBS -N ', "$job1\n";
    print $JOB1 '#PBS -q main', "\n";
    print $JOB1 '#PBS -M ems394@cornell.edu', "\n";
    print $JOB1 '#PBS -m abe', "\n";
    print $JOB1 '#PBS -A ged-intel11', "\n";
    print $JOB1 '#PBS -r n', "\n";
    print $JOB1 '#PBS -V', "\n";
    print $JOB1 "cd /mnt/home/emsch/work/Acey/ng_revision_jul2014/repeats ;\n";
    print $JOB1 "build_lmer_table -sequence $input_genome -freq $stem.freqs.out ;\n";
    print $JOB1 "qsub $job2 ;\n";

    close $JOB1;

    open my $JOB2, '>', $job2;

    print $JOB2 '#!/bin/bash -login', "\n";
    print $JOB2 '#PBS -l walltime=012:00:00', "\n";
    print $JOB2 '#PBS -l nodes=1:ppn=1', "\n";
    print $JOB2 '#PBS -l mem=50gb', "\n";
    print $JOB2 '#PBS -N ', "$job2\n";
    print $JOB2 '#PBS -q main', "\n";
    print $JOB2 '#PBS -M ems394@cornell.edu', "\n";
    print $JOB2 '#PBS -m abe', "\n";
    print $JOB2 '#PBS -A ged-intel11', "\n";
    print $JOB2 '#PBS -r n', "\n";
    print $JOB2 '#PBS -V', "\n";
    print $JOB2 "cd /mnt/home/emsch/work/Acey/ng_revision_jul2014/repeats ;\n";
    print $JOB2 "RepeatScout -sequence $input_genome -freq $stem.freqs.out -output $stem.reps.orig.fa ;\n";
    print $JOB2 "qsub $job3 ;\n";

    close $JOB2;

    open my $JOB3, '>', $job3;

    print $JOB3 '#!/bin/bash -login', "\n";
    print $JOB3 '#PBS -l walltime=001:00:00', "\n";
    print $JOB3 '#PBS -l nodes=1:ppn=1', "\n";
    print $JOB3 '#PBS -l mem=25gb', "\n";
    print $JOB3 '#PBS -N ', "$job3\n";
    print $JOB3 '#PBS -q main', "\n";
    print $JOB3 '#PBS -M ems394@cornell.edu', "\n";
    print $JOB3 '#PBS -m abe', "\n";
    print $JOB3 '#PBS -A ged-intel11', "\n";
    print $JOB3 '#PBS -r n', "\n";
    print $JOB3 '#PBS -V', "\n";
    print $JOB3 "cd /mnt/home/emsch/work/Acey/ng_revision_jul2014/repeats ;\n";
    print $JOB3 "filter-stage-1.prl $stem.reps.orig.fa > $stem.reps.fa ;\n";

    close $JOB3;
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

