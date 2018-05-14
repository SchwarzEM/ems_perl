#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;

my $dir = cwd;

my %ensembl = ( 'Homo_sapiens.GRCh37.70.pep.all.fa' => 1,
                'Mus_musculus.GRCm38.70.pep.all.fa' => 1, );

my @inputs = @ARGV;

foreach my $infile (@inputs) { 
    my $stem      = $infile;
    my $threshold = 1e-05;
    if ( $ensembl{$infile} ) { 
        $threshold = 1e-03;
    }
    $stem =~ s/\.fa\w*\z//;

    my $jobfile = 'job_' . $stem . '.hmmsearch.' . $threshold . '.pfam-A-26.0_19feb2013.sh';
    $jobfile    = safename($jobfile);

    open my $JOB, '>', $jobfile or die "Can't open job file $jobfile: $!";

    print $JOB '#!/bin/bash -login', "\n";
    print $JOB '#PBS -l walltime=012:00:00', "\n";
    print $JOB '#PBS -l nodes=1:ppn=8', "\n";
    print $JOB '#PBS -l mem=10gb', "\n";
    print $JOB '#PBS -N ', "$jobfile\n";
    print $JOB '#PBS -q main', "\n";
    print $JOB '#PBS -M ems394@cornell.edu', "\n";
    print $JOB '#PBS -m abe', "\n";
    print $JOB '#PBS -A ged-intel11', "\n";
    print $JOB '#PBS -r n', "\n";
    print $JOB '#PBS -V', "\n";
    print $JOB "cd $dir ;\n";
    print $JOB "module load HMMER/3.0 ;\n";
    print $JOB "hmmsearch --cpu 8 --noali -E $threshold -o /dev/null",
               " --tblout $stem.hmmsearch.$threshold.pfam-A-26.0.tblout.txt /mnt/scratch/emsch/dbs/Pfam-A.hmm $infile ;\n",
               ;

    close $JOB or die "Can't close filehandle to job file $jobfile: $!";
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

