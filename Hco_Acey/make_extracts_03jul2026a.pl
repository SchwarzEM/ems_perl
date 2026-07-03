#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;

my $targets = q{};
my $i       = 0;
my $batch   = q{};
my $build   = q{};

$targets = $ARGV[0] if $ARGV[0];
$i       = $ARGV[1] if $ARGV[1];

if ( (! $targets ) or (! $i ) ) {
    die "Format: make_extracts_03jul2026a.pl [target list] [start integer] => [script commands and script-building commands]\n";
}

open my $TARGETS, '<', $targets;
while ( my $target = <$TARGETS> ) {
    chomp $target;
    my $dir  = dirname($target);

    my $j  = sprintf "%02u", $i;
    $batch = "job_Necator_extract_2026.07.03.$j.sh";
    $build = "build_job_Necator_extract_2026.07.03.$j.sh";   

    open my $JOB1, '>', $batch;
    print $JOB1 '#!/bin/bash', "\n";
    print $JOB1 '#SBATCH --nodes=1', "\n";
    print $JOB1 '#SBATCH --partition=RM-shared', "\n";
    print $JOB1 '#SBATCH --time=04:00:00', "\n";
    print $JOB1 '#SBATCH --ntasks-per-node=2', "\n";
    print $JOB1 "#SBATCH --job-name=$batch\n";
    print $JOB1 '#SBATCH --mail-type=ALL', "\n";
    print $JOB1 "cd $dir ;\n";
    close $JOB1;

    open my $JOB2, '>', $build;
    print $JOB2 '$PROJECT/ems_perl/orthomcl_and_phylo/extract_syntelogs_26jun2026a.pl ';
    print $JOB2 "$target ";
    print $JOB2 '$PROJECT/necator/2026.01.14/genespace_04/annots/Aroian_genespace_04_seqids.2026.06.30.01.tsv.txt ';
    print $JOB2 '$PROJECT/necator/2026.01.14/genespace_04/annots/taxon2_pep.fasta_2026.06.30.01.tsv.txt ';
    print $JOB2 'pep.fa ';
    print $JOB2 '1>>';
    print $JOB2 "$batch ";
    print $JOB2 '2>';
    print $JOB2 "test$j.err ;";
    print $JOB2 "\n";
    close $JOB2;

    $i++;
    $j     = sprintf "%02u", $i;
    $batch = "job_Necator_extract_2026.07.03.$j.sh";
    $build = "build_job_Necator_extract_2026.07.03.$j.sh";

    open my $JOB3, '>', $batch;
    print $JOB3 '#!/bin/bash', "\n";
    print $JOB3 '#SBATCH --nodes=1', "\n";
    print $JOB3 '#SBATCH --partition=RM-shared', "\n";
    print $JOB3 '#SBATCH --time=04:00:00', "\n";
    print $JOB3 '#SBATCH --ntasks-per-node=2', "\n";
    print $JOB3 "#SBATCH --job-name=$batch\n";
    print $JOB3 '#SBATCH --mail-type=ALL', "\n";
    print $JOB3 "cd $dir ;\n";
    close $JOB3;

    open my $JOB4, '>', $build;
    print $JOB4 '$PROJECT/ems_perl/orthomcl_and_phylo/extract_syntelogs_26jun2026a.pl ';
    print $JOB4 "$target ";
    print $JOB4 '$PROJECT/necator/2026.01.14/genespace_04/annots/Aroian_genespace_04_seqids.2026.06.30.01.tsv.txt ';
    print $JOB4 '$PROJECT/necator/2026.01.14/genespace_04/annots/taxon2_dna.fasta_2026.06.30.01.tsv.txt ';
    print $JOB4 'dna.fa ';
    print $JOB4 '1>>';
    print $JOB4 "$batch ";
    print $JOB4 '2>';
    print $JOB4 "test$j.err ;";
    print $JOB4 "\n";
    close $JOB4;

    $i++;
}
close $TARGETS;


