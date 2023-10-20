#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;
use File::Basename;

my $start_dir = getcwd;
my $header    = 1;

while (my $reads = <>) {
    chomp $reads;
    if (! -r $reads) {
        die "Cannot access putative read file: $reads\n";
    }
    if ( $reads =~ /\.gz\z/xms ) {
        die "Should not try to compress an already compressed file: $reads\n";
    }
    my $stats = "$reads.seqkit_stats.txt";
    $stats    = basename($stats);

    if ($header) {
        print '#!/bin/bash', "\n";
        print '#SBATCH --nodes=1', "\n";
        print '#SBATCH --partition=RM-shared', "\n";
        print '#SBATCH --time=48:00:00', "\n";
        print '#SBATCH --ntasks-per-node=16', "\n";
        print '#SBATCH --job-name=job_Acey_seqkit.stats_2023.10.20.01.sh', "\n";
        print '#SBATCH --mail-type=ALL', "\n";
        print "cd $start_dir ;\n";
        print 'source $HOME/.bashrc_mamba ;', "\n";
        print '. $PROJECT/mambaforge-pypy3/etc/profile.d/mamba.sh ;', "\n";
        print 'mamba activate seqkit_2.4.0 ;', "\n";
        $header = 0;
    }
    print "seqkit stats --all --threads 16 $reads > $stats ;\n";
    print "gzip -1 $reads ;\n";
    print "ln -s $reads.gz ;\n";
}

