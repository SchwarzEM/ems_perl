#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;
use File::Basename;
use Scalar::Util qw(looks_like_number);

my $infile   = q{};
my $script   = q{};
my $curr_dir = getcwd;
my $date     = q{2001.01.01.01};
my $i        = 0;

$infile = $ARGV[0] if $ARGV[0];
$date   = $ARGV[1] if $ARGV[1];

# For this last variable, if I'm specifying it, I need to do two steps:
if ($ARGV[2]) {
    $i = $ARGV[2];
    $i--;
}

if ( (! $infile ) or ( $i < 0 ) or (! looks_like_number($i) ) or (! $i == int($i) ) ) {
    die "Format: make_bcf_bzip.scripts_09mar2026a.pl\n",
        "    [input table of uncompressed files]\n",
        "    [optional: date, default \"2001.01.01.01\"; should be clean date string]\n",
        "    [optional: starting index number for scripts; default 1, must be nonnegative integer]\n",
        "    => output batch scripts\n",
        ;
}

open my $INFILE, '<', $infile;
while ( my $input = <$INFILE> ) {
    chomp $input;
    my $output = q{};
    my $j      = q{};

    if ( $input =~ /\A \S+ \z/xms ) {
        $output = "$input.gz";
        if (-e $output ) {
            die "Intended script would overwrite existing file: $output\n";
        }
    }
    else {
        die "Cannot parse input $input\n";
    }

    $i++;
    $j = sprintf ("%03u", $i);

    $script = "job_bcf_bzip_" . "$date.$j.sh";
    $script = safename($script);

    open my $SCRIPT, '>', $script;

    print $SCRIPT '#!/bin/bash', "\n";
    print $SCRIPT '#SBATCH --nodes=1', "\n";
    print $SCRIPT '#SBATCH --partition=RM-shared', "\n";
    print $SCRIPT '#SBATCH --time=012:00:00', "\n";
    print $SCRIPT '#SBATCH --ntasks-per-node=4', "\n";
    print $SCRIPT "#SBATCH --job-name=$script\n";
    print $SCRIPT '#SBATCH --mail-type=ALL', "\n";
    print $SCRIPT "cd $curr_dir ;\n";
    print $SCRIPT 'source $HOME/.bashrc_mamba ;', "\n";
    print $SCRIPT '. $PROJECT/mambaforge-pypy3/etc/profile.d/mamba.sh ;', "\n";
    print $SCRIPT 'mamba activate samtools_1.17 ;', "\n";
    print $SCRIPT "bgzip -l 1 -@ 4 $input ;\n";
    print $SCRIPT "mamba deactivate ;\n";
    print $SCRIPT 'mamba activate bcftools_1.23 ;', "\n";
    print $SCRIPT "bcftools index --threads 4 $output ;\n";
    print $SCRIPT "mamba deactivate ;\n";

    close $SCRIPT;
}
close $INFILE;

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

