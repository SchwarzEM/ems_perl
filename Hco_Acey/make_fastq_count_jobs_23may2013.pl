#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use List::MoreUtils qw{ uniq };
use Cwd;

my $dir     = getcwd;
my @infiles = ();

while (my $input = <>) { 
    chomp $input;
    if (! -e $input) { 
        die "Can't find file: $input\n";
    }
    push @infiles, $input;
}

@infiles = uniq @infiles;

foreach my $infile (@infiles) { 
    my $basename = basename $infile;
    my $stem     = $basename;
    $stem        =~ s/\.fq\.gz\z//;
    my $job      = 'job_count_' . $stem . '_23may2013';
    my $job_file = $job . '.sh';
    $job_file    = safename($job_file);
    open my $JOB, '>', $job_file or die "Can't open job file $job_file: $!";

    print $JOB '#!/bin/bash -login', "\n";
    print $JOB '#PBS -l walltime=004:00:00', "\n";
    print $JOB '#PBS -l nodes=1:ppn=1', "\n";
    print $JOB '#PBS -l mem=10gb', "\n";
    print $JOB '#PBS -N ', "$job\n";
    print $JOB '#PBS -q main', "\n";
    print $JOB '#PBS -M ems394@cornell.edu', "\n";
    print $JOB '#PBS -m abe', "\n";
    print $JOB '#PBS -A ged-intel11', "\n";
    print $JOB '#PBS -r n', "\n";
    print $JOB '#PBS -V', "\n";
    print $JOB 'cd ', "$dir ;\n";
    print $JOB "zcat $infile | count_simple_fastq_residues.pl > $infile.count.txt ;\n";

    close $JOB or die "Can't close filehandle to job file $job_file: $!";
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

