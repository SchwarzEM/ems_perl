#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;
use File::Spec::Functions;  # catdir, catfile

my $infile_list = q{};
$infile_list = $ARGV[0] if $ARGV[0];

if (! $infile_list ) {
    die "Format: make_job_nippo_fastp_2023.09.14.01.pl [infile list] > [batch job of fastp against infiles]\n";
}

my $header      = 1;
my $start_dir = getcwd;

open my $LIST, '<', $infile_list;
while (my $infile_1 = <$LIST>) {
    chomp $infile_1;

    # Sample input from a file listing R1 readfiles:
    # /ocean/projects/mcb190015p/schwarze/nippo_steiner/nippo/genbank/uploads/ucr/SAMN37187589.IJ_02/IJ_rep01.R1.fastq.gz
    if ( $infile_1 =~ /\A (\S+) \/ ( [^\s\/]+ ) \/ ( [^\s\/]+ ) \.fastq \.gz \z/xms ) {
        my $dir1_stem = $1;
        my $dir2_stem = $2;
        my $file1_stem = $3;

        my $outdir = $dir2_stem;
        # remove any _02 suffixes on the input dir names
        $outdir    =~ s/_02\z//;

        if (! -d $outdir ) {
            die "Cannot see output directory: $outdir\n";
        }

        my $json    = "$file1_stem.filt1.json";
        my $html    = "$file1_stem.filt1.html";
        $json       = catfile($outdir, $json);
        $html       = catfile($outdir, $html);

        my $infile_2   = $infile_1;
        $infile_2      =~ s/\.R1\.fastq\.gz\z/.R2.fastq.gz/;

        my $file2_stem = $file1_stem;
        $file2_stem    =~ s/R1\z/R2/;

        if ( $infile_1 eq $infile_2 ) {
            die "Failure to distinguish $infile_1 from $infile_2\n";
        }
        if ( $file1_stem eq $file2_stem ) {
            die "Failure to distinguish $file1_stem from $file2_stem\n";
        }

        if (! -r $infile_1) {
            die "Cannot find or read infile_1: $infile_1\n";
        }
        if (! -r $infile_2) {
            die "Cannot find or read infile_2: $infile_2\n";
        }

        my $outfile_1  = "$file1_stem.filt1.fq";
        my $outfile_2  = "$file2_stem.filt1.fq";
        $outfile_1     = catfile($outdir, $outfile_1);
        $outfile_2     = catfile($outdir, $outfile_2);

        # print header lines only once, at the top of the output script

        print '#!/bin/bash', "\n" if $header;
        print '#SBATCH --nodes=1', "\n" if $header;
        print '#SBATCH --partition=RM-shared', "\n" if $header;
        print '#SBATCH --time=48:00:00', "\n" if $header;
        print '#SBATCH --ntasks-per-node=16', "\n" if $header;
        print '#SBATCH --job-name=job_nippo_fastp_2023.09.14.XX.sh', "\n" if $header;
        print '#SBATCH --mail-type=ALL', "\n" if $header;
        print "cd $start_dir ;\n" if $header;
        print 'source $HOME/.bashrc_mamba ;', "\n" if $header;
        print '. $PROJECT/mambaforge-pypy3/etc/profile.d/mamba.sh ;', "\n" if $header;
        print 'mamba activate fastp_0.23.2 ;', "\n" if $header;

        $header = 0;

        print "fastp --thread 16 --dont_overwrite --detect_adapter_for_pe --max_len1 50 --length_required 40 --json $json --html $html";
        print " --in1 $infile_1 --in2 $infile_2 --out1 $outfile_1 --out2 $outfile_2 ;\n";
    }

    else {
        die "Cannot parse input line: $infile_1\n";
    }
}
close $LIST;

print 'mamba deactivate ;', "\n" unless $header;
