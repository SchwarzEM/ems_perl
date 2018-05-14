#!/usr/bin/env perl

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    my $stem = $input;
    $stem    =~ s/\.fa\z//;

    my $target_dir = $stem . '.dir';
    my $qsub_script = 'job_iprscan_' . $stem . '_16nov2012.sh';
    system "mkdir $target_dir";
    system "mv $input $target_dir";
    open my $SCRIPT, '>', $qsub_script or die "Can't open qsub script: $qsub_script\n";

    print $SCRIPT '#!/bin/bash -login', "\n";
    print $SCRIPT '#PBS -l walltime=024:00:00', "\n";
    print $SCRIPT '#PBS -l nodes=1:ppn=4', "\n";
    print $SCRIPT '#PBS -l mem=10gb', "\n";
    print $SCRIPT "#PBS -N $qsub_script\n";
    print $SCRIPT '#PBS -q main', "\n";
    print $SCRIPT '#PBS -M ems394@cornell.edu', "\n";
    print $SCRIPT '#PBS -m abe', "\n";
    print $SCRIPT '#PBS -A ged-intel11', "\n";
    print $SCRIPT '#PBS -r n', "\n";
    print $SCRIPT '#PBS -V', "\n";
    print $SCRIPT "cd /mnt/scratch/emsch/blast2go/iprscan/$target_dir ;\n";
    print $SCRIPT 'module load iprscan ;', "\n";
    print $SCRIPT 'iprscan -cli -email ems394@cornell.edu', 
                  " -i $input -o $stem.iprscan.raw -format raw -iprlookup -goterms ;\n", 
                  ;

    close $SCRIPT or die "Can't close filehandle to qsub script: $qsub_script\n";
    system "mv $qsub_script $target_dir";
} 

