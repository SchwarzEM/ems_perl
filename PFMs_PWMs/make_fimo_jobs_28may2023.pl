#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;
use File::Basename;
use File::Spec::Functions;  # catdir, catfile

my $start_dir = getcwd;

my $list     = q{};
my $motif_db = q{};
my $bgfile   = q{};
my $gen_dna  = q{};
my $script   = q{};
my $header   = 1;

$list     = $ARGV[0] if $ARGV[0];
$motif_db = $ARGV[1] if $ARGV[1];
$bgfile   = $ARGV[2] if $ARGV[2];
$gen_dna  = $ARGV[3] if $ARGV[3];
$script   = $ARGV[4] if $ARGV[4];

if ( (! $list ) or (! $motif_db ) or (! $bgfile ) or (! $gen_dna ) or (! $script ) ) {
    die "Format: make_fimo_jobs_27may2023.pl",
        " [search motif list] [single motif file] [Markov model file] [DNA FASTA file] [script name]",
        "\n",
        ;
}

if (! -r $motif_db ) {
    die "Cannot read putative motif database file: $motif_db\n";
}
if (! -r $bgfile) {
    die "Cannot read putative Markov model background file: $bgfile\n";
}
if (! -r $gen_dna) {
    die "Cannot read putative genomic DNA file: $gen_dna\n";
}

open my $LIST, '<', $list;
while (my $motif = <$LIST>) {
    chomp $motif;
    if ( $motif !~ /\A \S+ \z/xms ) {
        die "Motif in list $list has a forbidden space character: \"$motif\"\n";
    }

    my $fimo_dir  = $motif . '_fimo';
    $fimo_dir     = safename($fimo_dir);

    # Only print this once at the start of the output.
    if ($header) {
        print '#!/bin/bash', "\n";
        print '#SBATCH --nodes=1', "\n";
        print '#SBATCH --partition=RM-shared', "\n";
        print '#SBATCH --time=048:00:00', "\n";
        print '#SBATCH --ntasks-per-node=4', "\n";
        print '#SBATCH --job-name=', "$script\n";
        print '#SBATCH --mail-type=ALL', "\n";
        print "cd $start_dir ;\n";
        print 'source $HOME/.bashrc_mamba ;', "\n";
        print '. $PROJECT/mambaforge-pypy3/etc/profile.d/mamba.sh ;', "\n";
        print "mamba activate meme_5.4.1 ;\n";
        $header = 0;
    }

    print "fimo --motif $motif --bgfile $bgfile ";
    print "--max-stored-scores 10000000 -o $fimo_dir --thresh 1e-4 $motif_db ";
    print "$gen_dna ;\n";
}
close $LIST;

if (! $header) {
    print "mamba deactivate ;\n";
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
