#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;
use File::Basename;
use Scalar::Util qw(looks_like_number);

my $infile = q{};
my $chroms = q{};
my $date   = q{2001.01.01.01};
my $i      = 0;

$infile = $ARGV[0] if $ARGV[0];
$chroms = $ARGV[1] if $ARGV[1];
$date   = $ARGV[2] if $ARGV[2];

# For this last variable, if I'm specifying it, I need to do two steps:
if ($ARGV[3]) {
    $i = $ARGV[3];
    $i--;
}

if ( (! $infile ) or (! $chroms ) or ( $i < 0 ) or (! looks_like_number($i) ) or (! $i == int($i) ) ) {
    die "Format: make_freebayes.scripts_04mar2026a.pl\n",
        "    [input table of BAMS and their target genomes]\n",
        "    [list of individual chromosomes/scaffolds for which to generate VCFs]\n",
        "    [optional: date, default \"2001.01.01.01\"; should be clean date string]\n",
        "    [optional: index number, default 0, must be nonnegative integer]\n",
        "    => output batch scripts\n",
        ;
}

my $stem     = q{};
my $script   = q{};
my $j        = 0;
my $in_bam   = q{};
my $genome   = q{};

my $curr_dir = getcwd;

my @chrs = ();  # previously was hard-coded as: qw( Necator_chrI Necator_chrII Necator_chrIII Necator_chrIV Necator_chrV Necator_chrX );

open my $CHROMS, '<', $chroms;
while ( my $input = <$CHROMS> ) {
    chomp $input;
    if ( $input =~ /\A \S+ \z/xms ) {
        push @chrs, $input;
    }
    else {
        die "From target chromosome/scaffold list file $chroms, cannot parse: $input\n";
    }
}
close $CHROMS;

my $chrom_count = @chrs;
$chrom_count--;

open my $INFILE, '<', $infile;
while ( my $input = <$INFILE> ) {
    chomp $input;
    if ( $input =~ /\A ((\S+) \.bam) \t (\S+) \z/xms ) {
        $in_bam = $1;
        $stem   = $2;
        $genome = $3;
    }
    else {
        die "Cannot parse input $input\n";
    }

    foreach my $k (0..$chrom_count) {
        $i++;
        $j = sprintf ("%02u", $i);
        my $chr     = $chrs[$k];
        my $out_vcf = "$stem.$chr.fb.$date.vcf";
        $out_vcf    = basename($out_vcf);
        $out_vcf    = safename($out_vcf);

        $script = "job_Necator_freebayes.etc_" . "$date.$j.sh";
        $script = safename($script);
        open my $SCRIPT, '>', $script;

        print $SCRIPT '#!/bin/bash', "\n";
        print $SCRIPT '#SBATCH --nodes=1', "\n";
        print $SCRIPT '#SBATCH --partition=RM-shared', "\n";
        print $SCRIPT '#SBATCH --time=048:00:00', "\n";
        print $SCRIPT '#SBATCH --ntasks-per-node=4', "\n";
        print $SCRIPT "#SBATCH --job-name=$script\n";
        print $SCRIPT '#SBATCH --mail-type=ALL', "\n";
        print $SCRIPT "cd $curr_dir ;\n";
        print $SCRIPT 'source $HOME/.bashrc_mamba ;', "\n";
        print $SCRIPT '. $PROJECT/mambaforge-pypy3/etc/profile.d/mamba.sh ;', "\n";

        print $SCRIPT "mamba activate freebayes_1.3.10 ;\n";
        print $SCRIPT "freebayes -f $genome -r $chr --standard-filters --report-monomorphic -g 200 $in_bam > $out_vcf ;\n";
        print $SCRIPT 'mamba deactivate ;', "\n";

        close $SCRIPT;
    }
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

