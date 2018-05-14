#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $i = 0;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \.func_bool.txt \z/xms ) { 
        my $stem = $1;
        my $output = "$stem.func_bool.analysis";
        $output    = safename($output);

        $i++;
        my $j = sprintf "%02i", $i;

        my $job = 'job_func_hyper_' . $stem . "_2016.02.06.$j.sh";
        $job    = safename($job);

        open my $JOB, '>', $job;

        print $JOB '#!/bin/bash -login', "\n";
        print $JOB '#PBS -l walltime=001:00:00', "\n";
        print $JOB '#PBS -l nodes=1:ppn=1', "\n";
        print $JOB '#PBS -l mem=16gb', "\n";
        print $JOB "#PBS -N $job\n";
        print $JOB '#PBS -q main', "\n";
        print $JOB '#PBS -M ems394@cornell.edu', "\n";
        print $JOB '#PBS -m abe', "\n";
        print $JOB '#PBS -A ged', "\n";
        print $JOB '#PBS -r n', "\n";
        print $JOB '#PBS -V', "\n";
        print $JOB 'cd /mnt/home/emsch/work/2015/adrienne/go_analysis ;', "\n";
        print $JOB "mkdir $output ;\n";
        print $JOB "func_hyper -i $input -t go_201512-termdb-tables -o $output ;\n";

        close $JOB;
    }
    else {
        die "Cannot parse input: $input\n";
    }
}


sub safename {
    my $_filename = $_[0];
    my $_orig_filename = $_filename;
    if (-e $_orig_filename) {
        my $_suffix1 = 1;
        $_filename = $_filename . ".$_suffix1";
        while (-e $_filename) {
            $_suffix1++;
            $_filename =~ s/\.\d+\z//xms;
            $_filename = $_filename . ".$_suffix1";
        }
    }
    return $_filename;
}

