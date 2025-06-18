#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $in_sam = q{};
my $stem   = q{};
my $script = q{};
my $i      = 0;
my $j      = 0;

while ( $in_sam = <> ) {
    chomp $in_sam;
    if ( $in_sam =~ /\A (\S+) \.sam \z/xms ) {
        $stem = $1;
    }
    else {
        die "Cannot parse infile $in_sam\n";
    }
    my $out_bam1 = "$stem.orig.bam";
    my $out_bam2 = "$stem.sort.bam";
    my $out_bam3 = "$stem.dup_flag.bam";
    my $out_bam4 = "$stem.dup_flag.sort.bam";
    my $tmp_dir  = 'tmp_' . $stem;

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


    print $SCRIPT "mamba activate samtools_1.17 ;\n";
    print $SCRIPT "samtools view -@ 15 -b $in_sam > $out_bam1 ;\n";
    print $SCRIPT "samtools sort -@ 15 -o $out_bam2 $out_bam1 ;\n";
    print $SCRIPT "samtools index -@ 15 $out_bam2 ;\n";
    print $SCRIPT 'mamba deactivate ;', "\n";

    print $SCRIPT "mamba activate sambamba_1.0.1 ;\n";
    print $SCRIPT "sambamba markdup -t 16 --tmpdir=", "$tmp_dir $out_bam2 $out_bam3 ;\n";
    print $SCRIPT 'mamba deactivate ;', "\n";

    print $SCRIPT "mamba activate samtools_1.17 ;\n";
    print $SCRIPT "samtools sort -@ 15 -o $out_bam4 $out_bam3 ;\n";
    print $SCRIPT "samtools index -@ 15 $out_bam4 ;\n";
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

