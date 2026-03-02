#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Getopt::Long;
use Cwd;
use List::MoreUtils qw(uniq);
use File::Basename;

my $infile      = q{};
my $init_dir    = getcwd;
my $i           = 1;

my $help;

GetOptions ( 'infiles=s' => \$infile,
             'start=i'   => \$i,
             'help'      => \$help,   );

if ( $help or (! $infile) ) {
    die "Format: make_necator_samtools_01mar2026.pl\n",
        "    --infile|-i   <list of BAM files to sort and index; with full filenames>\n",
        "    --start|-s    <starting index for script numbers; default of 1>\n",
        "    --help|-h     [print this message]\n",
        ;
}

$i--;

open my $INFILE, '<', $infile;

while (my $input = <$INFILE>) {
    chomp $input;

    if (! -e $input ) {
        die "In input file list $infile, cannot find: $input\n";
    }
    my $basename = basename($input);
    my $work_dir = dirname($input);
    my $output   = q{};

    if ( $basename =~ /\A (\S+) \. bam \z/xms ) {
        my $stem = $1;
        $output  = "$stem.sort.bam";
    }
    else {
        die "Cannot parse basename in: $input\n";
    }

    $i++;
    my $j = sprintf "%03u", $i;

    my $output_script = "job_necator_samtools_2026.03.02.$j.sh";
    $output_script    = safename($output_script);

    open my $OUT, '>', $output_script;

    print $OUT '#!/bin/bash', "\n";
    print $OUT '#SBATCH --nodes=1', "\n";
    print $OUT '#SBATCH --partition=RM-shared', "\n";
    print $OUT '#SBATCH --time=24:00:00', "\n";
    print $OUT '#SBATCH --ntasks-per-node=16', "\n";
    print $OUT '#SBATCH --job-name=', $output_script, "\n";
    print $OUT '#SBATCH --mail-type=ALL', "\n";

    print $OUT "cd $init_dir ;\n";
    print $OUT "cd $work_dir ;\n";
    print $OUT 'source $HOME/.bashrc_mamba ;', "\n";
    print $OUT '. $PROJECT/mambaforge-pypy3/etc/profile.d/mamba.sh ;', "\n";
    print $OUT 'mamba activate samtools_1.17 ;', "\n";
    print $OUT "samtools sort -@ 15 -o $output $basename ;\n";
    print $OUT "samtools index -@ 15 $output ;\n";
    print $OUT "mamba deactivate ;\n";

    close $OUT;
}

sub safename {
    my $_filename = $_[0];
    my $_orig_filename = $_filename;
    if (-e $_orig_filename) {
        my $_suffix1 = 1;
        $_filename = $_filename . ".$_suffix1";
        while (-e $_filename) {
            $_suffix1++;
            $_filename =~ s/\.\d+\z//xms;
            $_filename = $_filename . ".$_suffix1";
        }
    }
    return $_filename;
}

