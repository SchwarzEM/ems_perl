#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

my $umass = q{};
my $washu = q{};

$umass = $ARGV[0] if $ARGV[0];
$washu = $ARGV[1] if $ARGV[1];

if ( (! $umass ) or (! $washu ) ) {
    die "Format: comp_motif_qvals_10oct2025.pl [UMass motif qvals] [WashU motif qvals] > [UMass/WashU motif qval comparison]\n";
}

my $header = "Motif\tUMass_qval\tWashU-only_qval\tWashU-only.to.UMessRatio";
my $data_ref;

open my $UMASS, '<', $umass;
while ( my $input = <$UMASS> ) {
    chomp $input;
    if ( $input =~ /\A ([^\t]+) \t .* \t ([^\t]+) \z/xms ) { 
        my $motif = $1;
        my $qval  = $2;
        if ( (! looks_like_number($qval) ) and ( $qval ne 'q-value' ) ) {
            die "From data file $umass, cannot parse q-value as number: $qval\n";
        }
        $data_ref->{'motif'}->{$motif}->{umass_qval} = $qval;
    }
    else {
        die "From data file $umass, cannot parse: $input\n";
    }
}
close $UMASS;

open my $WASHU, '<', $washu;
while ( my $input = <$WASHU> ) {
    chomp $input;
    if ( $input =~ /\A ([^\t]+) \t .* \t ([^\t]+) \z/xms ) {
        my $motif = $1;
        my $qval  = $2;
        if ( (! looks_like_number($qval) ) and ( $qval ne 'q-value' ) ) {
            die "From data file $umass, cannot parse q-value as number: $qval\n";
        }
        $data_ref->{'motif'}->{$motif}->{washu_qval} = $qval;
    }
    else {
        die "From data file $umass, cannot parse: $input\n";
    }
}
close $WASHU;

my @raw_motifs = keys %{ $data_ref->{'motif'} };
foreach my $raw_motif (@raw_motifs) {
    my $umass_qval = 1;
    my $washu_qval = 1;
    if ( ( exists $data_ref->{'motif'}->{$raw_motif}->{umass_qval} ) and ( $data_ref->{'motif'}->{$raw_motif}->{umass_qval} ne 'q-value' ) )  {
        $umass_qval = $data_ref->{'motif'}->{$raw_motif}->{umass_qval};
    }
    if ( ( exists $data_ref->{'motif'}->{$raw_motif}->{washu_qval} ) and ( $data_ref->{'motif'}->{$raw_motif}->{washu_qval} ne 'q-value' ) ) {
        $washu_qval = $data_ref->{'motif'}->{$raw_motif}->{washu_qval};
    }
    my $ratio = ($washu_qval / $umass_qval);

    print "$header\n" if $header;
    $header = q{};

    if ( ( $umass_qval < 0.1 ) or ( $washu_qval < 0.1 ) ) {
        print "$raw_motif\t$umass_qval\t$washu_qval\t$ratio\n";
    }
}

