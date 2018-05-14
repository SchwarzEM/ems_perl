#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use File::Basename;

my $curr_dir = cwd;

my @fastqs = ();

while (my $input = <>) {
    chomp $input;
    push @fastqs, $input;
}

foreach my $fastq (@fastqs) { 
    my $basename = basename $fastq;
    my $job_name = $basename ;

    $job_name     =~ s/\.gz\z//;
    my $stem_name = $job_name ;

    $job_name =~ s/\.fq\z//;
    $job_name =~ s/\.fastq\z//;
    $job_name = 'job_count_' . $job_name . '.sh';
    $job_name = safename($job_name);

    open my $JOB, '>', $job_name or die "Can't open job script $job_name\n";

    my $header =   '#!/bin/bash -login'         . "\n"
                 . '#PBS -l walltime=012:00:00' . "\n"
                 . '#PBS -l nodes=1:ppn=1'      . "\n"
                 . '#PBS -l mem=60gb'           . "\n"
                 . '#PBS -N '                   . "$job_name \n"
                 . '#PBS -q main'               . "\n"
                 . '#PBS -M ems394@cornell.edu' . "\n"
                 . '#PBS -m abe'                . "\n"
                 . '#PBS -A ged-intel11'        . "\n"
                 . '#PBS -r n' . "\n"
                 . '#PBS -V'   . "\n"
                 . "cd $curr_dir ;\n"
                 ;

    print $JOB $header;
    print $JOB "zcat $fastq | fastq2fa_simple.pl | count_fasta_residues.pl -i - -e > modEnc_count_", $stem_name, "_13feb2013.txt ;\n";
    print $JOB "\n";
    close $JOB or die "Can't close filehandle to job script $job_name\n";
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

