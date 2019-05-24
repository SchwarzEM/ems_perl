#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;
use File::Spec::Functions;  # catdir, catfile

my $header      = 1;
my $start_dir = getcwd;

while (my $input = <>) {
    chomp $input;

    # Sample input directory name:
    # /pylon5/mc5phvp/schwarze/Acey/2019.05/Acey_RNAseq_data/filt1_RNAseq_reads/17751_Acey_intestines_19dpi_noDEX_biorep_1.filt1.fq
    if ( $input =~ /\A (\S+) \/ ([^\s\/]+) \z/xms ) {
        my $dir    = $1;
        my $subdir = $2;
        my $outdir = catdir($dir, 'filt2_RNAseq_reads');
        if (! -r $outdir ) {
            die "Cannot read out: $outdir\n";
        }
        my $outfile = "$subdir.filt2.fq.gz";
        $outfile    = catfile($outdir, $outfile);

        # print header lines only once, at the top of the output script

        print '#!/bin/bash', "\n" if $header;
        print '#SBATCH --nodes=1', "\n" if $header;
        print '#SBATCH --partition=RM-shared', "\n" if $header;
        print '#SBATCH --time=48:00:00', "\n" if $header;
        print '#SBATCH --ntasks-per-node=1', "\n" if $header;
        print '#SBATCH --constraint=EGRESS', "\n" if $header;
        print '#SBATCH --job-name=job_acey_fastp_2019.05.24.02.sh', "\n" if $header;
        print '#SBATCH --mail-type=ALL', "\n" if $header;
        print "cd $start_dir ;\n" if $header;
        print '. $SCRATCH/anaconda2/etc/profile.d/conda.sh ;', "\n" if $header;
        print 'conda activate fastp_0.20.0 ;', "\n" if $header;

        $header = 0;

        print "zcat $input", '/*.gz', ' | fastp --dont_overwrite';
        print ' --adapter_fasta $SCRATCH/Acey/2019.05/Acey_RNAseq_data/indiv_seqs/illumina_adaptors_all_2016.10.08.01.all.fa';
        print ' --n_base_limit 0 --max_len1 50';
        print " --stdin --out1 $outfile ;\n";
    }

    else {
	die "Cannot parse input line: $input\n";
    }
}

print 'conda deactivate ;', "\n" unless $header;
