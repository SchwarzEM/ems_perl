#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $pgs  = q{};
my $alns = q{};

$pgs  = $ARGV[0] if $ARGV[0];
$alns = $ARGV[1] if $ARGV[1];

if ( (! $pgs ) or (! $alns ) ) {
    die "Format: pick_pgene_aligns_21jul2026a.pl [selected pangene table] [alignment table] > [selected alignment table]\n";
}

my %pick = ();

open my $PGS, '<', $pgs;
while ( my $input = <$PGS> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t \S+ \z/xms ) {
        my $pg = $1;
        $pick{$pg} = 1;
    }
    else {
        die "From pangene table $pgs, cannot parse: $input\n";
    }
}
close $PGS;

open my $ALNS, '<', $alns;
while ( my $input = <$ALNS> ) {
    chomp $input;
    # Sample input:
    # syntelogs_01/pangene_1001.dna.mafft_aln.trimal_aut1.fa
    if ( $input =~ /\A \S+ (pangene_\d+) \. \S+ \z/xms ) {
        my $pg = $1;
        if ( exists $pick{$pg} ) {
            print "$input\n";
        }
    }
    else {
        die "From alignment list $alns, cannot parse: $input\n";
    }
}
close $ALNS;
