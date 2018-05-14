#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use File::Basename;

my @input_files = @ARGV;
my $motif_table = shift @input_files;
my $q_val_prog  = shift @input_files;

foreach my $input_file (@input_files) {
    # Sample:  rsem_log10_slices/Acey_2012.10.24.pgene_log10.24.PI.vs.L3i_17apr2014.txt
    my $basename = basename($input_file);
    if ( $basename =~ /\A Acey_2012\.10\.24\.pgene_log10\.(\S+)_17apr2014\.txt \z/xms ) { 
        my $stem = $1;

        my $qsub_script = 'job_wilcox_motifs_log.' . $stem . '_' . '20apr2014.sh';
        open my $SCRIPT, '>', $qsub_script;

        print $SCRIPT '#!/bin/bash -login', "\n";
        print $SCRIPT '#PBS -l walltime=012:00:00', "\n";
        print $SCRIPT '#PBS -l nodes=1:ppn=1', "\n";
        print $SCRIPT '#PBS -l mem=10gb', "\n";
        print $SCRIPT "#PBS -N $qsub_script\n";
        print $SCRIPT '#PBS -q main', "\n";
        print $SCRIPT '#PBS -M ems394@cornell.edu', "\n";
        print $SCRIPT '#PBS -m abe', "\n";
        print $SCRIPT '#PBS -A ged-intel11', "\n";
        print $SCRIPT '#PBS -r n', "\n";
        print $SCRIPT '#PBS -V', "\n";  
        print $SCRIPT "cd /mnt/home/emsch/work/Acey/ng_revision_mar2014/rev_table_s10 ;\n";
        print $SCRIPT "/mnt/home/emsch/perl.svn/trunk/Hco_Acey/wilcoxon_gene_annots_w_q_vals.pl -a $motif_table -r $input_file -q $q_val_prog",
                      ' 1>', "Acey_2012.10.24.motifs_wilcoxon_log10.$stem", "_20apr2014.txt",
                      ' 2>', "Acey_2012.10.24.motifs_wilcoxon_log10.$stem", "_20apr2014.err",
                      " ;\n",
                      ;

        close $SCRIPT;
    }
    else { 
        die "Can't parse input file: $input_file\n";
    }
} 

