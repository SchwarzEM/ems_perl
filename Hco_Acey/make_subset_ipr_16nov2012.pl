#!/usr/bin/env perl

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (\S+\.)(\d+) \z/xms ) { 
    my $stem  = $1;
    my $digit = $2;
    my $index = $digit;
    $index    =~ s/\A[0]//;

    my $index_next       = $index;
    $index_next++;
    my $digit_next       = sprintf "%02i", $index_next;
    my $qsub_script      = 'job_iprscan_' . $stem . $digit      . '_16nov2012.sh';
    my $qsub_script_next = 'job_iprscan_' . $stem . $digit_next . '_16nov2012.sh';

    open my $SCRIPT, '>', $qsub_script or die "Can't open qsub script: $qsub_script\n";

    print $SCRIPT '#!/bin/bash -login', "\n";
    print $SCRIPT '#PBS -l walltime=024:00:00', "\n";
    print $SCRIPT '#PBS -l nodes=1:ppn=4', "\n";
    print $SCRIPT '#PBS -l mem=20gb', "\n";
    print $SCRIPT "#PBS -N $qsub_script\n";
    print $SCRIPT '#PBS -q main', "\n";
    print $SCRIPT '#PBS -M ems394@cornell.edu', "\n";
    print $SCRIPT '#PBS -m abe', "\n";
    print $SCRIPT '#PBS -A ged-intel11', "\n";
    print $SCRIPT '#PBS -r n', "\n";
    print $SCRIPT '#PBS -V', "\n";
    print $SCRIPT "cd /mnt/scratch/emsch/blast2go/iprscan ;\n";
    print $SCRIPT 'module load iprscan ;', "\n";
    print $SCRIPT 'iprscan -cli -email ems394@cornell.edu', 
                  " -i $input -o $input.iprscan.raw -format raw -iprlookup -goterms ;\n", 
                  ;
    print $SCRIPT "qsub $qsub_script_next ;\n";

    close $SCRIPT or die "Can't close filehandle to qsub script: $qsub_script\n";
    }
    else { 
        die "Can't parse input: $input\n";
    }
} 

