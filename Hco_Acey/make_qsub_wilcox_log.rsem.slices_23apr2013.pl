#!/usr/bin/env perl

use strict;
use warnings;

my @input_files = @ARGV;
my $motif_table = shift @input_files;

foreach my $input_file (@input_files) {
    if ( $input_file =~ /\A Acey_2012\.10\.24\.pgene_log10\.(\S+)_24mar2013\.txt \z/xms ) { 
        my $stem = $1;

        my $qsub_script = 'job_wilcox_motifs_log.' . $stem . '_' . '23apr2013.sh';
        open my $SCRIPT, '>', $qsub_script or die "Can't open qsub script: $qsub_script\n";

        print $SCRIPT '#!/bin/bash -login', "\n";
        print $SCRIPT '#PBS -l walltime=008:00:00', "\n";
        print $SCRIPT '#PBS -l nodes=1:ppn=1', "\n";
        print $SCRIPT '#PBS -l mem=10gb', "\n";
        print $SCRIPT "#PBS -N $qsub_script\n";
        print $SCRIPT '#PBS -q main', "\n";
        print $SCRIPT '#PBS -M ems394@cornell.edu', "\n";
        print $SCRIPT '#PBS -m abe', "\n";
        print $SCRIPT '#PBS -A ged-intel11', "\n";
        print $SCRIPT '#PBS -r n', "\n";
        print $SCRIPT '#PBS -V', "\n";  
        print $SCRIPT "cd /mnt/lustre_scratch_2012/emsch/Acey.3/wilcoxon_motifs ;\n";
        print $SCRIPT "wilcoxon_gene_annots.pl -a $motif_table -r $input_file",
                      ' 1>', "Acey_2012.10.24.motifs_wilcoxon_log10.$stem", "_23apr2013.txt",
                      ' 2>', "Acey_2012.10.24.motifs_wilcoxon_log10.$stem", "_23apr2013.err",
                      " ;\n",
                      ;

        close $SCRIPT or die "Can't close filehandle to qsub script: $qsub_script\n";
    }
    else { 
        die "Can't parse input file: $input_file\n";
    }
} 

