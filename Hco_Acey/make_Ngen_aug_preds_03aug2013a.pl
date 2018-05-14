#!/usr/bin/env perl

use strict;
use warnings;

my @genome_stems = qw ( brenneri_WS235.HM.2013.05.05
                        c_angaria.WS235.genomic
                        Csp11_WUv2.dna
                        C_sp1_scaffolds_17jul2013.tag
                        Csp5_WS235.HM.2013.05.05
                        japonica_WS235.HM.2013.05.05
                        remanei_WS235.HM.2013.05.05 );

# Omitted from the list, this time: C_sp16_scaffolds_17jul2013.tag and Csp9_Berkeley_30nov2012_k65.GC.HMerg.tag

foreach my $gen_stem (@genome_stems) { 
    my $script_stem = 'job_augustus_' . $gen_stem . '.03aug2013a' ;
    my $script      = $script_stem . '.sh';
    $script         = safename($script);

    open my $SCRIPT, '>', $script or die "Can't open script $script: $!";

    print $SCRIPT '#!/bin/bash -login', "\n", ;
    print $SCRIPT '#PBS -l walltime=072:00:00', "\n", ;
    print $SCRIPT '#PBS -l nodes=1:ppn=1', "\n", ;
    print $SCRIPT '#PBS -l mem=10gb', "\n", ;
    print $SCRIPT '#PBS -N ', "$script_stem\n";
    print $SCRIPT '#PBS -q main', "\n", ;
    print $SCRIPT '#PBS -M ems394@cornell.edu', "\n", ;
    print $SCRIPT '#PBS -m abe', "\n", ;
    print $SCRIPT '#PBS -A ged-intel11', "\n", ;
    print $SCRIPT '#PBS -r n', "\n", ;
    print $SCRIPT '#PBS -V', "\n", ;
    print $SCRIPT 'cd /mnt/scratch/emsch/Ngen ;', "\n", ;
    print $SCRIPT 'augustus --strand=both --genemodel=partial --singlestrand=false',
                  ' --alternatives-from-sampling=true --minexonintronprob=0.1 --minmeanexonintronprob=0.4',
                  ' --uniqueGeneId=true --protein=on --cds=on --UTR=on --introns=on --start=on --stop=on --cds=on --codingseq=on',
                  ' --species=caenorhabditis --extrinsicCfgFile=/mnt/home/emsch/src/augustus.2.7/config/extrinsic/extrinsic.ME.cfg',
                  ' --progress=true',
                  ' --outfile=', "$gen_stem.03aug2013.aug $gen_stem.fa ;\n",
                  ;

    close $SCRIPT or die "Can't close filehandle to script $script: $!";
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

