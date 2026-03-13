#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

my $infile = q{};
my $bin    = q{};

$infile = $ARGV[0] if $ARGV[0];
$bin    = $ARGV[1] if $ARGV[1];

if ( (! $infile ) or (! $bin ) or (! looks_like_number($bin) ) ) {
    die "Format: snpden2bed_13mar2026.pl [input snpden file] [bin size] > [output BED file]\n";
}

open my $INFILE, '<', $infile;

while ( my $input = <$INFILE> ) {
    chomp $input;

    # Sample inputs:
    # CHROM       BIN_START       SNP_COUNT       VARIANTS/KB
    # Ilik2_chrI      0       892     8.92
    # Ilik2_chrI      100000  1088    10.88
    # Ilik2_chrI      200000  1384    13.84

    if ( $input =~ /\A (\S+) \t (\S+) \t \S+ \t (\S+) \z/xms ) {
        my $chrom   = $1;
        my $coord   = $2;
        my $density = $3;
        if ( ( $chrom ne 'CHROM' ) and ( looks_like_number($coord) ) and ( looks_like_number($density) ) ) {
            my $start = $coord;
            my $end   = ($coord + $bin);
            print "$chrom\t$start\t$end\t$density\n";
        }
    }
    else {
        die "Can't parse input line: $input\n";
    }
}

close $INFILE;

