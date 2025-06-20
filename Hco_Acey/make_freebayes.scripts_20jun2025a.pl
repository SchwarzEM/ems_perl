#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $stem    = q{};
my $script  = q{};
my $i       = 0;
my $j       = 0;
my $in_bam  = q{};
my $genome  = q{};

while ( my $input = <> ) {
    chomp $input;
    if ( $input =~ /\A ((\S+) \.bam) \t (\S+) \z/xms ) {
        $in_bam = $1;
        $stem   = $2;
        $genome = $3;
    }
    else {
        die "Cannot parse input $input\n";
    }

    $i++;
    $j = sprintf ("%02u", $i);

    $script = "job_Necator_freebayes.etc_2025.06.20.$j.sh";
    $script = safename($script);
    open my $SCRIPT, '>', $script;

    print $SCRIPT '#!/bin/bash', "\n";
    print $SCRIPT '#SBATCH --nodes=1', "\n";
    print $SCRIPT '#SBATCH --partition=RM-shared', "\n";
    print $SCRIPT '#SBATCH --time=024:00:00', "\n";
    print $SCRIPT '#SBATCH --ntasks-per-node=16', "\n";
    print $SCRIPT "#SBATCH --job-name=$script\n";
    print $SCRIPT '#SBATCH --mail-type=ALL', "\n";
    print $SCRIPT 'cd $PROJECT/necator/2023.09.12/cov/illu ;', "\n";
    print $SCRIPT 'source $HOME/.bashrc_mamba ;', "\n";
    print $SCRIPT '. $PROJECT/mambaforge-pypy3/etc/profile.d/mamba.sh ;', "\n";

    print $SCRIPT "mamba activate freebayes_1.3.10 ;\n";
    print $SCRIPT "freebayes --report-monomorphic -f $genome --standard-filters $in_bam > $stem.fb01.vcf ;\n";
    print $SCRIPT 'mamba deactivate ;', "\n";

    close $SCRIPT;
}

sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}

