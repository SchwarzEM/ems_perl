#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;
use File::Spec::Functions;  # catdir, catfile

my $header      = 1;
my $start_dir = getcwd;

while (my $infile_1 = <>) {
    chomp $infile_1;

    # Sample input from 'ls $PROJECT/Acey/2019.05/pub_RNAseq_data_03/filt0_RNAseq_reads/*_1.fq';
    # /ocean/projects/mcb190015p/schwarze/Acey/2019.05/pub_RNAseq_data_03/filt0_RNAseq_reads/WashU_adult_fem_01_1.fastq

    if ( $infile_1 =~ /\A (\S+) \/ filt0_RNAseq_reads \/ (\S+) _1\.fq \z/xms ) {
        my $dir_stem  = $1;
        my $file_stem = $2;

        my $outdir = catdir($dir_stem, 'filt2_all.nt_RNAseq_reads');
        if (! -r $outdir ) {
            die "Cannot read out: $outdir\n";
        }

        my $infile_2   = $infile_1;
        $infile_2      =~ s/_1\.fq\z/_2.fq/;

       	if (! -r $infile_1) {
       	    die	"Cannot find or read infile_1: $infile_1\n";
        }
        if (! -r $infile_2) {
            die "Cannot find or read infile_2: $infile_2\n";
        }

	my $outfile_1  = "$file_stem.filt2.all.nt_1.fq";
        my $outcount_1 = "$file_stem.filt2.all.nt_1.fq.readcount.txt";

        my $outfile_2  = "$file_stem.filt2.all.nt_2.fq";
        my $outcount_2 = "$file_stem.filt2.all.nt_2.fq.readcount.txt";

        my $json       = "$file_stem.filt2.all.nt.json";
        my $html       = "$file_stem.filt2.all.nt.html";

        $outfile_1     = catfile($outdir, $outfile_1);
        $outcount_1    = catfile($outdir, $outcount_1);

        $outfile_2     = catfile($outdir, $outfile_2);
        $outcount_2    = catfile($outdir, $outcount_2);

        # print header lines only once, at the top of the output script

        print '#!/bin/bash', "\n" if $header;
        print '#SBATCH --nodes=1', "\n" if $header;
        print '#SBATCH --partition=RM-shared', "\n" if $header;
        print '#SBATCH --time=48:00:00', "\n" if $header;
        print '#SBATCH --ntasks-per-node=24', "\n" if $header;
        print '#SBATCH --job-name=job_acey_fastp_20XX.XX.XX.XX.sh', "\n" if $header;
        print '#SBATCH --mail-type=ALL', "\n" if $header;
        print "cd $start_dir ;\n" if $header;
        print '. $PROJECT/anaconda3/etc/profile.d/conda.sh ;', "\n" if $header;
        print 'conda activate fastp_0.20.0 ;', "\n" if $header;

        $header = 0;

        print "fastp --thread 24 --dont_overwrite --json $json --html $html";
        print ' --n_base_limit 0 --length_required 50';
        print " --in1 $infile_1 --in2 $infile_2 --out1 $outfile_1 --out2 $outfile_2 ;\n";
        print "cat $outfile_1 | count_simple_fastq_residues.pl > $outcount_1 ;\n";
        print "cat $outfile_2 | count_simple_fastq_residues.pl > $outcount_2 ;\n";
    }

    else {
	die "Cannot parse input line: $infile_1\n";
    }
}

print 'conda deactivate ;', "\n" unless $header;
