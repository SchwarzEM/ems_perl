#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $stem    = q{};
my $script  = q{};
my $i       = 22;
my $j       = 0;
my $genome  = q{};

while ( my $input = <> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \.sam \t (\S+) \z/xms ) {
        $stem   = $1;
        $genome = $2;
    }
    else {
        die "Cannot parse input $input\n";
    }
    my $in_bam1 = "$stem.sort.bam";
    my $in_bam2 = "$stem.dup_flag.sort.bam";

    $i++;
    $j = sprintf ("%02u", $i);

    $script = "job_Necator_freebayes.etc_2025.06.18.$j.sh";
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
    print $SCRIPT "freebayes -f $genome -C 5 --standard-filters $in_bam2 > $stem.dup_flag.sort.fb01.vcf ;\n";
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

