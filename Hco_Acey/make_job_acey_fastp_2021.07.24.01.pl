#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;
use File::Spec::Functions;  # catdir, catfile

my $header      = 1;
my $start_dir = getcwd;

while (my $infile = <>) {
    chomp $infile;

    # Sample input from 'ls $PROJECT/Acey/2019.05/Acey_RNAseq_data/filt1_RNAseq_reads/*.filt1.fq':
    # /ocean/projects/mcb190015p/schwarze/Acey/2019.05/Acey_RNAseq_data/filt1_RNAseq_reads/17751_Acey_intestines_19dpi_noDEX_biorep_1.filt1.fq
    if ( $infile =~ /\A (\S+) \/ filt1_RNAseq_reads \/ (\S+) \.filt1\.fq \z/xms ) {
        my $dir_stem  = $1;
        my $file_stem = $2;
        my $outdir = catdir($dir_stem, 'filt2_RNAseq_reads');
        if (! -r $outdir ) {
            die "Cannot read out: $outdir\n";
        }
        my $outfile = "$file_stem.filt2.fq";
        my $json    = "$file_stem.filt2.json";
        my $html    = "$file_stem.filt2.html";
        $outfile    = catfile($outdir, $outfile);

        # print header lines only once, at the top of the output script

        print '#!/bin/bash', "\n" if $header;
        print '#SBATCH --nodes=1', "\n" if $header;
        print '#SBATCH --partition=RM-shared', "\n" if $header;
        print '#SBATCH --time=48:00:00', "\n" if $header;
        print '#SBATCH --ntasks-per-node=16', "\n" if $header;
        print '#SBATCH --job-name=job_acey_fastp_2021.07.24.01.sh', "\n" if $header;
        print '#SBATCH --mail-type=ALL', "\n" if $header;
        print "cd $start_dir ;\n" if $header;
        print '. $PROJECT/anaconda3/etc/profile.d/conda.sh ;', "\n" if $header;
        print 'conda activate fastp_0.20.0 ;', "\n" if $header;

        $header = 0;

        print "fastp --thread 16 --dont_overwrite --detect_adapter_for_pe --json $json --html $html";
        print ' --n_base_limit 0 --length_required 50';
        print " --in1 $infile --out1 $outfile ;\n";
    }

    else {
	die "Cannot parse input line: $infile\n";
    }
}

print 'conda deactivate ;', "\n" unless $header;
