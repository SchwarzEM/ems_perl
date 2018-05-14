#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $qsub_start_script = 'run_bash5tools_jobs_02sep2014.sh';
$qsub_start_script    = safename($qsub_start_script);

open my $QSUB, '>', $qsub_start_script;

my $i = 0;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A ( \/mnt\/home\/emsch\/work\/Ngen\/pacbio\/meyer\/ [^\s\/]+ \/ [^\s\/]+ \/ Analysis_Results ) \/ ( ( [^\s\/]+ ) \.bas\.h5 ) \z/xms ) {
        my $working_directory = $1;
        my $input_bas_h5_file = $2;
        my $prefix            = $3;

        $i++;
        my $j = sprintf "%02i", $i;
        my $job_name = "job_bash5tools_02sep2014.$j";
        my $job_file = "$working_directory/$job_name.sh";
        $job_file = safename($job_file);

        open my $JOB, '>', $job_file;

        print $JOB '#!/bin/bash -login', "\n";
        print $JOB '#PBS -l walltime=024:00:00', "\n";
        print $JOB '#PBS -l nodes=1:ppn=1', "\n";
        print $JOB '#PBS -l mem=25gb', "\n";
        print $JOB "#PBS -N $job_name\n";
        print $JOB '#PBS -q main', "\n";
        print $JOB '#PBS -M ems@emstech.org', "\n";
        print $JOB '#PBS -m abe', "\n";
        print $JOB '#PBS -A ged-intel11', "\n";
        print $JOB '#PBS -r n', "\n";
        print $JOB '#PBS -V', "\n";
        print $JOB "cd $working_directory ;\n";
        print $JOB 'module load Python/2.7.2 ;', "\n";
        print $JOB 'module load h5py/2.1.3 ;', "\n";
        print $JOB 'module load NumPy/1.8.1 ;', "\n";
        print $JOB "bash5tools.py --outFilePrefix $prefix.subreads --readType subreads --outType fastq --minReadScore 0.75 --minLength 100 $input_bas_h5_file ;\n";

        close $JOB;

        print $QSUB "    qsub $job_file ;\n";
    }
    else { 
        die "This highly specific, hard-coded Perl jobbuilder script cannot parse input line: $input\n";
    }
}

close $QSUB;

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

