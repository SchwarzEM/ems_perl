#!/usr/bin/env perl

use strict;
use warnings;

my @k_vals = (33..48);

foreach my $k_val (@k_vals) { 
    my $script = 'job_27feb2013_Acey_and_Cel_pme_TPM.5plus_reads.k' . $k_val . '.27feb2013a.sh' ;
    open my $SCRIPT, '>', $script or die "Can't open script $script: $!";

    print $SCRIPT '#!/bin/bash -login', "\n", ;
    print $SCRIPT '#PBS -l walltime=018:00:00', "\n", ;
    print $SCRIPT '#PBS -l nodes=1:ppn=1', "\n", ;
    print $SCRIPT '#PBS -l mem=10gb', "\n", ;
    print $SCRIPT '#PBS -N job_27feb2013_Acey_and_Cel_pme_TPM.5plus_reads.k', $k_val, '.27feb2013a.sh', "\n", ;
    print $SCRIPT '#PBS -q main', "\n", ;
    print $SCRIPT '#PBS -M ems394@cornell.edu', "\n", ;
    print $SCRIPT '#PBS -m abe', "\n", ;
    print $SCRIPT '#PBS -A ged-intel11', "\n", ;
    print $SCRIPT '#PBS -r n', "\n", ;
    print $SCRIPT '#PBS -V', "\n", ;
    print $SCRIPT 'cd /mnt/lustre_scratch_2012/emsch/rsem_work/2013.02.23 ;', "\n", ;
    print $SCRIPT 'module load cluster3/1.50 ;', "\n", ;
    print $SCRIPT 'cluster -f Acey_and_Cel_modencode_pme_TPM.5plus_reads.27feb2013.strict_orths.pzero_TPMs.txt -l -g 7 -k ', 
                  $k_val, 
                  ' -r 10000 -u Acey_and_Cel.27feb2013a.k', 
                  "$k_val ;\n", 
                  ;

    close $SCRIPT or die "Can't close filehandle to script $script: $!";
}

