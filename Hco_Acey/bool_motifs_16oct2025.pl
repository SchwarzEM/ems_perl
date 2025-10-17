#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $infile = q{};
my $motif  = q{};

$infile = $ARGV[0] if $ARGV[0];
$motif  = $ARGV[1] if $ARGV[1];

if ( (! $infile ) or (! $motif ) ) {
    die "Format: bool_motifs_16oct2025.pl [gene-motif TSV] [single_motif] > [Booleanized gene-single_motif TSV]\n";
}

open my $INFILE, '<', $infile;
while ( my $input = <$INFILE> ) {
    my $output = q{};
    chomp $input;
    if ( $input =~ /\A Gene \t [^\t]+ \z/xms ) {
        $output = "Gene\t$motif";
    }
    elsif ( $input =~ /\A (\S+) \t [^\t]* $motif [^\t]* \z/xms ) {
        my $gene = $1;
        $output = "$gene\ttrue";
    }
    elsif ( $input =~ /\A (\S+) \t [^\t]* \z/xms ) {
        my $gene = $1;
        $output = "$gene\tfalse";
    }
    else {
        die "From infile $infile, cannot parse: $input\n";
    }
    print "$output\n";
}
close $INFILE;



