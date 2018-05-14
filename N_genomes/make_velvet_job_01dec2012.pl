#!/usr/bin/env perl

use strict;
use warnings;
use List::MoreUtils qw(uniq);
use Scalar::Util qw(looks_like_number);

my @initial_k_vals = @ARGV;
my @k_vals         = sort { $a <=> $b } @initial_k_vals;
@k_vals            = uniq @k_vals;

my $val_count = @k_vals;
if ( $val_count == 0 ) { 
    die "Format: make_velvet_job_01dec2012.pl [k value or values]\n";
}

foreach my $k_val (@k_vals) { 
    # Reality-check the numbers.
    if (! looks_like_number($k_val) ) { 
        die "k_val $k_val does not look like number\n";
    }
    if ( $k_val != int($k_val) ) { 
        die "k_val $k_val does not look like integer\n";
    }
    if ( $k_val <= 2 ) { 
        die "k_val $k_val is less than 3\n";
    }

    my $velveth = 'job_velveth_Csp9_Berkeley_01dec2012_k' . $k_val . '.sh';
    my $velvetg = 'job_velvetg_Csp9_Berkeley_01dec2012_k' . $k_val . '.sh';

    open my $VELVETH, '>', $velveth or die "Can't open velveth script $velveth\n";

    print $VELVETH '#!/bin/bash -login', "\n", ;
    print $VELVETH '#PBS -l walltime=096:00:00', "\n", ;
    print $VELVETH '#PBS -l nodes=1:ppn=16', "\n", ;
    print $VELVETH '#PBS -l mem=200gb', "\n", ;
    print $VELVETH '#PBS -N job_velveth_Csp9_Berkeley_30nov2012_k', "$k_val \n", ;
    print $VELVETH '#PBS -q main', "\n", ;
    print $VELVETH '#PBS -M ems394@cornell.edu', "\n", ;
    print $VELVETH '#PBS -m abe', "\n", ;
    print $VELVETH '#PBS -A ged-intel11', "\n", ;
    print $VELVETH '#PBS -r n', "\n", ;
    print $VELVETH '#PBS -V', "\n", ;
    print $VELVETH 'cd /mnt/scratch/emsch/sp9 ;', "\n", ;
    print $VELVETH 'velveth Csp9_Berkeley_30nov2012_k',
                   "$k_val $k_val",
                   ' -create_binary -fasta.gz',
                   ' -short CS004_NoIndex_L002_R12.khmer142-1.NONcont.se.fa.gz',
                   ' -shortPaired CS004_NoIndex_L002_R12.khmer142-1.NONcont.pe.fa.gz ;',
                   "\n",
                   ;
    print $VELVETH "qsub $velvetg ;\n", ;

    close $VELVETH or die "Can't close filehandle to velveth script $velveth\n";

    open my $VELVETG, '>', $velvetg or die "Can't open velvetg script $velvetg\n";

    print $VELVETG '#!/bin/bash -login ', "\n", ;
    print $VELVETG '#PBS -l walltime=096:00:00', "\n", ;
    print $VELVETG '#PBS -l nodes=1:ppn=16', "\n", ;
    print $VELVETG '#PBS -l mem=200gb', "\n", ;
    print $VELVETG '#PBS -N job_velvetg_Csp9_Berkeley_30nov2012_k', "$k_val\n", ;
    print $VELVETG '#PBS -q main', "\n", ;
    print $VELVETG '#PBS -M ems394@cornell.edu', "\n", ;
    print $VELVETG '#PBS -m abe', "\n", ;
    print $VELVETG '#PBS -A ged-intel11', "\n", ;
    print $VELVETG '#PBS -r n', "\n", ;
    print $VELVETG '#PBS -V', "\n", ;
    print $VELVETG 'cd /mnt/scratch/emsch/sp9 ;', "\n", ;
    print $VELVETG 'velvetg Csp9_Berkeley_30nov2012_k', "$k_val -cov_cutoff auto -exp_cov auto -min_contig_lgth 200 -ins_length 335 ;\n", ;

    close $VELVETG or die "Can't close filehandle to velvetg script $velvetg\n";
}

