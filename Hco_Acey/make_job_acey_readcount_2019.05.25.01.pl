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
    my $outfile = "$input.readcount.txt";
    
    print '#!/bin/bash', "\n" if $header;
    print '#SBATCH --nodes=1', "\n" if $header;
    print '#SBATCH --partition=RM-shared', "\n" if $header;
    print '#SBATCH --time=48:00:00', "\n" if $header;
    print '#SBATCH --ntasks-per-node=1', "\n" if $header;
    print '#SBATCH --constraint=EGRESS', "\n" if $header;
    print '#SBATCH --job-name=job_acey_readcount_2019.05.25.01.sh', "\n" if $header;
    print '#SBATCH --mail-type=ALL', "\n" if $header;
    print "cd $start_dir ;\n" if $header;

    $header = 0;

    print "zcat $input | count_simple_fastq_residues.pl > $outfile ;\n";
}
