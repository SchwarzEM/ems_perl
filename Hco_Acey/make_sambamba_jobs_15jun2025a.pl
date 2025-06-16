#!/usr/bin/env perl 

use strict;
use warnings;
use autodie;

my $infile = q{};
my $xx     = q{};
my $stem   = q{};

$infile = $ARGV[0] if $ARGV[0];
$xx     = $ARGV[1] if $ARGV[1];

if ( $infile =~ /(\S+) \.bam/xms ) {
    $stem = $1;
}

my $tmpdir = 'tmp_' . $stem . '_dir';

if ( $infile and $stem and $xx ) {
    print '#!/bin/bash', "\n";
    print '#SBATCH --nodes=1', "\n";
    print '#SBATCH --partition=RM-shared', "\n";  
    print '#SBATCH --time=001:00:00', "\n";
    print '#SBATCH --ntasks-per-node=16', "\n";
    print "#SBATCH --job-name=job_Necator_sambamba_2025.06.15.$xx.sh\n";
    print '#SBATCH --mail-type=ALL', "\n";
    print 'cd $PROJECT/necator/2023.09.12/cov/illu ;', "\n";
    print 'source $HOME/.bashrc_mamba ;', "\n";
    print '. $PROJECT/mambaforge-pypy3/etc/profile.d/mamba.sh ;', "\n";
    print 'mamba activate sambamba_1.0.1 ;', "\n";
    print "sambamba markdup -t 16 --tmpdir=$tmpdir $infile $stem.pre.bam ;\n";
    print 'mamba deactivate ;', "\n";
    print 'mamba activate samtools_1.17 ;', "\n";
    print "samtools sort -@ 15 -o $stem.flag01.bam $stem.pre.bam ;\n";
    print "samtools index -@ 15 $stem.flag01.bam ;\n";
    print 'mamba deactivate ;', "\n";
    print 'mamba activate sambamba_1.0.1 ;', "\n";
    print "sambamba depth window -w 1000000 -t 16 --min-coverage=0 -o $stem.flag01.depth.txt $stem.flag01.bam ;\n";
    print "sambamba flagstat -b -t 16 $stem.flag01.bam > $stem.flag01.flagstat.01.txt ;\n";
    print 'mamba deactivate ;', "\n";
}

