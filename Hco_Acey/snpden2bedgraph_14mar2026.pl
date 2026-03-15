#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

my $infile = q{};
my $bin    = q{};
my $chr    = q{};

$infile = $ARGV[0] if $ARGV[0];
$bin    = $ARGV[1] if $ARGV[1];
$chr    = $ARGV[2] if $ARGV[2];

my %chr_max = ();

if ( (! $infile ) or (! $bin ) or (! looks_like_number($bin) ) or (! $chr ) ) {
    die "Format: snpden2bedgraph_14mar2026.pl [input snpden file] [bin size] [chr size file] > [output BEDGraph file]\n";
}

open my $CHR, '<', $chr;
while ( my $input = <$CHR> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\d+) \z/xms ) {
        my $chr_id   = $1;
        my $chr_size = $2;
        if (! looks_like_number($chr_size) ) {
            die "From chromosome size file $chr, cannot identify chromosome nt size in: $input\n";
        }
        $chr_max{$chr_id} = $chr_size;
    }
    else {
        die "From chromosome size file $chr, cannot parse: $input\n";
    }
}
close $CHR;

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
            my $start  = $coord;
            my $end    = ($coord + $bin);
            if (! exists  $chr_max{$chrom} ) {
                die "Can't identify nt size for chromosome $chrom\n";
            }
            my $max_nt = $chr_max{$chrom};
            if ( $end > $max_nt ) {
                $end = $max_nt;
            }
            print "$chrom\t$start\t$end\t$density\n";
        }
    }
    else {
        die "Can't parse input line: $input\n";
    }
}
close $INFILE;

