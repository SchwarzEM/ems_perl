#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

my $input2job =	q{};
my $genemap   = q{};
my $salmon_db = q{};
my $boots     = q{};

$input2job    = $ARGV[0] if $ARGV[0];
$genemap      = $ARGV[1] if $ARGV[1];
$salmon_db    = $ARGV[2] if $ARGV[2];
$boots        = $ARGV[3] if $ARGV[3];

if (! $input2job ) {
    die "Format: make_salmon_jobs_06sep2020_v02.pl [input2job file] [cds2gene file] [salmon_db dir] [bootstrap integer] > [job.sh]\n";
}

if (! -r $genemap) {
    die "Cannot read gene map file: $genemap\n";
}
if (! -d $salmon_db) {
    die "Putative salmon db is not an actual directory: $salmon_db\n";
}

if (! $boots ) {
    $boots = 0;
}
elsif ( (! looks_like_number($boots) ) or ( $boots < 0 ) or ( $boots != int($boots) ) ) {
    die "Bootstrap number must be a non-negative integer (0 is allowed): $boots\n";
}

open my $INPUT2JOB, '<', $input2job;

while (my $input = <$INPUT2JOB>) {
    chomp $input;
    if ( $input =~ /\A ((\S+)_1\.filt1\.fq\.gz) \t (\S+) \z/xms ) {
        my $infile1 = $1;
        my $stem    = $2;
        my $label   = $3;

        if (! -r $infile1) {
            die "Cannot read input file 1: $infile1\n";
        }

        my $infile2 = $stem . '_2.filt1.fq.gz';
        if (! -r $infile2) {
            die "Cannot read input file 2: $infile2\n";
        }

	my $outfile = $label . '_' . "$boots-boots.salmon";
        if (-e $outfile) {
            die "Planned output file already exists: $outfile\n";
        }         

        my $salmon_command = "salmon --no-version-check quant --threads 8 --libType A "
                             . "--seqBias --gcBias --posBias --validateMappings --rangeFactorizationBins 4 --numBootstraps $boots "
                             . "--geneMap $genemap "
                             . "--index $salmon_db "
                             . "--mates1 $infile1 "
                             . "--mates2 $infile2 "
                             . "--output $outfile ;"
                             ;

        print "$salmon_command\n";
    }
    else {
        die "Cannot parse input2job line: $input\n";
    }
}

close $INPUT2JOB;
