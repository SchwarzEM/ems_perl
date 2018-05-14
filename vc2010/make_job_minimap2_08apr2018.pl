#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = 
    '#!/bin/bash -login'. "\n".
    '#PBS -l walltime=024:00:00'. "\n".
    '#PBS -l nodes=1:ppn=8'. "\n".
    '#PBS -l mem=32gb'. "\n".
    '#PBS -N job_mass_n2_cds.to.vc2010_minimap2_2018.04.08.01.sh'. "\n".
    '#PBS -q main'. "\n".
    '#PBS -M ems394@cornell.edu'. "\n".
    '#PBS -m abe'. "\n".
    '#PBS -A ged'. "\n".
    '#PBS -r n'. "\n".
    '#PBS -V'. "\n" ;

my $workdir = '/mnt/home/emsch/work/VC2010/final_feb2018/minimap2/cds_dnas';
my $target  = '/mnt/home/emsch/work/VC2010/final_feb2018/seqs/vc2010.draft-20180305.pilon.fasta';

while (my $query = <>) {
    chomp $query;
    if ( $query =~ /\A (\S+)\.fa \z/xms ) { 
        my $stem = $1;
        my $output = "$stem.minimap2.vc2010.sam";
        my $error  = "$stem.minimap2.vc2010.err";

        if ($header) {
            print $header;
            print "cd $workdir ;\n";
            $header = q{};
        }

        # minimap2 -a -t 8 -x asm5  MtDNA.fa  chrM_pilon.fa    1>vc2010_M.minimap2.n2_MtDNA.sam  2>vc2010_M.minimap2.err ;
        print "minimap2 -a -t 8 -x asm5 $target $query 1>", "$output 2>", "$error ;\n";
    }
}

