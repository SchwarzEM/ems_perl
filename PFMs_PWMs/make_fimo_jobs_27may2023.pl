#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;
use File::Basename;
use File::Spec::Functions;  # catdir, catfile

my $start_dir = getcwd;

my $list    = q{};
my $bgfile  = q{};
my $gen_dna = q{};
my $date    = q{};

$list    = $ARGV[0] if $ARGV[0];
$bgfile  = $ARGV[1] if $ARGV[1];
$gen_dna = $ARGV[2] if $ARGV[2];
$date    = $ARGV[3] if $ARGV[3];

if ( (! $list ) or (! $bgfile ) or (! $gen_dna ) or (! $date ) ) {
    die "Format: make_fimo_jobs_27may2023.pl [motif list] [Markov model file] [DNA FASTA file] [date or other id]\n";
}

if (! -r $bgfile) {
    die "Cannot read putative Markov model background file: $bgfile\n";
}
if (! -r $gen_dna) {
    die "Cannot read putative genomic DNA file: $gen_dna\n";
}

open my $LIST, '<', $list;
while (my $input = <$LIST>) {
    chomp $input;

    # Enforce reality of input file:
    if (! -r $input) {
        die "From input list file $list, cannot read putative file: $input\n";
    }
    # Enforce likely usability of the target MEME file:
    if ( $input !~ / \.meme\.txt \z/xms ) { 
        die "From input list file $list, input file lacks suffix .meme.txt (and may not work for fimo): $input\n";
    }

    # Do *not* safename $work_dir because it is supposed to already exist, and we are going to go back into it!
    my $work_dir = dirname($input);
    my $stem     = basename($input);
    $stem        =~ s/\.meme\.txt//;

    my $script   = 'job_fimo_' . $stem . '_' . "$date.sh";
    $script      = safename($script);

    my $full_script = catfile($work_dir, $script);

    my $fimo_dir  = $stem . '_fimo_' . $date;
    $fimo_dir     = safename($fimo_dir);

    # strictly speaking unnecessary, but it makes it easier to keep things straight, so do this:
    my $meme_file = $input;

    open my $SCRIPT, '>', $script;

    print $SCRIPT '#!/bin/bash', "\n";
    print $SCRIPT '#SBATCH --nodes=1', "\n";
    print $SCRIPT '#SBATCH --partition=RM-shared', "\n";
    print $SCRIPT '#SBATCH --time=001:00:00', "\n";
    print $SCRIPT '#SBATCH --ntasks-per-node=1', "\n";
    print $SCRIPT '#SBATCH --job-name=', "$script\n";
    print $SCRIPT '#SBATCH --mail-type=ALL', "\n";

    print $SCRIPT "cd $start_dir ;\n";
    print $SCRIPT 'source $HOME/.bashrc_mamba ;', "\n";
    print $SCRIPT '. $PROJECT/mambaforge-pypy3/etc/profile.d/mamba.sh ;', "\n";
    print $SCRIPT "mamba activate meme_5.4.1 ;\n";

    print $SCRIPT "fimo --bgfile $bgfile ";
    print $SCRIPT "--max-stored-scores 10000000 -o $fimo_dir --thresh 1e-4 $meme_file ";
    print $SCRIPT "$gen_dna ;\n";

    print $SCRIPT "mamba deactivate ;\n";

    close $SCRIPT;
}

close $LIST;

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


