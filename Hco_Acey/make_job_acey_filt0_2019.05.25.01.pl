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

    if ( $input =~ /\A (\S+) \/ ([^\s\/]+) \z/xms ) {
        my $dir    = $1;
        my $subdir = $2;
        my $outdir = catdir($dir, 'filt0_RNAseq_reads');
        if (! -r $outdir ) {
            die "Cannot read out: $outdir\n";
        }
        my $outfile = "$subdir.filt0.fq";
        $outfile    = catfile($outdir, $outfile);

        # print header lines only once, at the top of the output script

        print '#!/bin/bash', "\n" if $header;
        print '#SBATCH --nodes=1', "\n" if $header;
        print '#SBATCH --partition=RM-shared', "\n" if $header;
        print '#SBATCH --time=48:00:00', "\n" if $header;
        print '#SBATCH --ntasks-per-node=1', "\n" if $header;
        print '#SBATCH --constraint=EGRESS', "\n" if $header;
        print '#SBATCH --job-name=job_acey_filt0_2019.05.25.01.sh', "\n" if $header;
        print '#SBATCH --mail-type=ALL', "\n" if $header;
        print "cd $start_dir ;\n" if $header;

        $header = 0;

        print "zcat $input", '/*.gz', " > $outfile ;\n";
    }

    else {
	die "Cannot parse input line: $input\n";
    }
}

