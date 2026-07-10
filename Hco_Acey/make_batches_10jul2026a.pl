#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;

my $max  = 0;
my $date = q{};
my $dir  = getcwd;

$max    = $ARGV[0] if $ARGV[0];
$date   = $ARGV[1] if $ARGV[1];

if ( (! $max ) or (! $date ) ) {
    die "Format: make_batches_10jul2026a.pl [N, maximum batch number] [date in script names] => [N starter batch job files]\n";
}

foreach my $i ( 1..$max ) {
    my $j = sprintf "%02u", $i;
    my $outfile = "job_Necator_mass.hyphy_$date.$j.sh";
    $outfile    = safename($outfile);

    open my $OUTFILE, '>', $outfile;

    print $OUTFILE '#!/bin/bash', "\n";
    print $OUTFILE '#SBATCH --nodes=1', "\n";
    print $OUTFILE '#SBATCH --partition=RM-shared', "\n";
    print $OUTFILE '#SBATCH --time=48:00:00', "\n";
    print $OUTFILE '#SBATCH --ntasks-per-node=16', "\n";
    print $OUTFILE "#SBATCH --job-name=$outfile\n";
    print $OUTFILE '#SBATCH --mail-type=ALL', "\n";
    print $OUTFILE "cd $dir ;\n";
    print $OUTFILE 'source $HOME/.bashrc_mamba ;', "\n";
    print $OUTFILE '. $PROJECT/mambaforge-pypy3/etc/profile.d/mamba.sh ;', "\n";
    print $OUTFILE 'mamba activate hyphy_2.5.100 ;', "\n";

    close $OUTFILE;
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

