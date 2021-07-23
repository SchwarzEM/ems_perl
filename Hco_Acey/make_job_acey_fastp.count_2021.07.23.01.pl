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

    # Sample input from 'ls $PROJECT/Acey/2019.05/pub_RNAseq_data_03/filt0_RNAseq_reads/*_1.fq';
    # /ocean/projects/mcb190015p/schwarze/Acey/2019.05/pub_RNAseq_data_03/filt0_RNAseq_reads/WashU_adult_fem_01_1.fastq

    if ( $infile =~ /\A (\S+) \/ filt0_RNAseq_reads \/ (\S+) _1\.fq \z/xms ) {
        my $dir_stem  = $1;
        my $file_stem = $2;
        my $outdir = catdir($dir_stem, 'filt2_50nt_RNAseq_reads');
        if (! -r $outdir ) {
            die "Cannot read out: $outdir\n";
        }
        my $outfile = "$file_stem.filt2.fq.gz";
        my $outcount = "$file_stem.filt2.fq.readcount.txt";
        my $json    = "$file_stem.filt2.json";
        my $html    = "$file_stem.filt2.html";
        $outfile    = catfile($outdir, $outfile);
        $outcount    = catfile($outdir, $outcount);

        # print header lines only once, at the top of the output script

        print '#!/bin/bash', "\n" if $header;
        print '#SBATCH --nodes=1', "\n" if $header;
        print '#SBATCH --partition=RM-shared', "\n" if $header;
        print '#SBATCH --time=48:00:00', "\n" if $header;
        print '#SBATCH --ntasks-per-node=8', "\n" if $header;
        print '#SBATCH --job-name=job_acey_fastp.count_20XX.XX.XX.XX.sh', "\n" if $header;
        print '#SBATCH --mail-type=ALL', "\n" if $header;
        print "cd $start_dir ;\n" if $header;
        print '. $PROJECT/anaconda3/etc/profile.d/conda.sh ;', "\n" if $header;
        print 'conda activate fastp_0.20.0 ;', "\n" if $header;

        $header = 0;

        print "fastp --thread 8 --dont_overwrite --json $json --html $html ";
        print ' --n_base_limit 0 --length_required 50 --max_len1 50';
        print " --in1 $infile --out1 $outfile ;\n";
        print "zcat $outfile | count_simple_fastq_residues.pl > $outcount ;\n";
    }

    else {
	die "Cannot parse input line: $infile\n";
    }
}

print 'conda deactivate ;', "\n" unless $header;
