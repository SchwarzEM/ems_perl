#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $stem    = q{};
my $script  = q{};
my $i       = 48;
my $j       = 0;
my $genome  = '/ocean/projects/mcb190015p/schwarze/necator/2022.07.29/assemblies/Necator_2022.05.29.02.chrs_only.fa';

my @chrs = qw( Necator_chrI Necator_chrII Necator_chrIII Necator_chrIV Necator_chrV Necator_chrX );

foreach my $k (0..5) {
        $i++;
        $j = sprintf ("%02u", $i);
        my $chr     = $chrs[$k];
        my $out_vcf = "Necator_2022.05.29.02.chrs.ngm.8_genos.$chr.fb01.vcf";
        $out_vcf    = safename($out_vcf);

        $script = "job_Necator_freebayes.etc_2025.06.21.$j.sh";
        $script = safename($script);
        open my $SCRIPT, '>', $script;

        print $SCRIPT '#!/bin/bash', "\n";
        print $SCRIPT '#SBATCH --nodes=1', "\n";
        print $SCRIPT '#SBATCH --partition=RM-shared', "\n";
        print $SCRIPT '#SBATCH --time=048:00:00', "\n";
        print $SCRIPT '#SBATCH --ntasks-per-node=8', "\n";
        print $SCRIPT "#SBATCH --job-name=$script\n";
        print $SCRIPT '#SBATCH --mail-type=ALL', "\n";
        print $SCRIPT 'cd $PROJECT/necator/2023.09.12/cov/illu ;', "\n";
        print $SCRIPT 'source $HOME/.bashrc_mamba ;', "\n";
        print $SCRIPT '. $PROJECT/mambaforge-pypy3/etc/profile.d/mamba.sh ;', "\n";

        print $SCRIPT "mamba activate freebayes_1.3.10 ;\n";
        print $SCRIPT "freebayes -f $genome -r $chr --standard-filters --report-monomorphic -g 200";
        print $SCRIPT ' Necator_2022.05.29.02.chrs.ngm.Aroian_illu.v_sens.dup_flag.sort.bam';
        print $SCRIPT ' Necator_2022.05.29.02.chrs.ngm.Anhui_illu.v_sens.dup_flag.sort.bam';
        print $SCRIPT ' Necator_2022.05.29.02.chrs.ngm.Ilik2_illu.v_sens.dup_flag.sort.bam';
        print $SCRIPT ' Necator_2022.05.29.02.chrs.ngm.Mag3_illu.v_sens.dup_flag.sort.bam';
        print $SCRIPT ' Necator_2022.05.29.02.chrs.ngm.Oita_illu.v_sens.dup_flag.sort.bam';
        print $SCRIPT ' Necator_2022.05.29.02.chrs.ngm.Swissman_illu.v_sens.dup_flag.sort.bam';
        print $SCRIPT ' Necator_2022.05.29.02.chrs.ngm.Ilik4_illu.sens.dup_flag.sort.bam';
        print $SCRIPT ' Necator_2022.05.29.02.chrs.ngm.obscurus_illu.sens.dup_flag.sort.bam';
        print $SCRIPT " > $out_vcf ;\n";
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

