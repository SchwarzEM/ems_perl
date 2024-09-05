#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $equivs = q{};
my $annots = q{};

$equivs = $ARGV[0] if $ARGV[0];
$annots = $ARGV[1] if $ARGV[1];

if ( (! $equivs) or (! $annots) ) {
    die "Format: st_genes2equivs_05sep2024.pl [equivalance groups] [gene-indexed annots] > [equiv-indexed annots]\n";
}

my $data_ref;

open my $EQUIVS, '<', $equivs;
while ( my $equiv = <$EQUIVS> ) {
    chomp $equiv;
    if ( $equiv =~ /\A \S .* \S \z/xms ) {
        my @genes = split '; ', $equiv;
        foreach my $gene (@genes) {
            if ( exists $data_ref->{'gene'}->{$gene}->{'equiv'} ) {
                die "Redundant mapping of gene $gene to equivalence group $data_ref->{'gene'}->{$gene}->{'equiv'}\n";
            }
            $data_ref->{'gene'}->{$gene}->{'equiv'} = $equiv;
        }
    }
    else {
        die "From equivalence group list $equivs, cannot parse: $equiv\n";
    }
}
close $EQUIVS;

open my $ANNOTS, '<', $annots;
while ( my $annot = <$ANNOTS> ) {
    chomp $annot;
    if ( $annot =~ /\A (\S+) \t (.+) \z/xms ) {
        my $gene = $1;
        my $data = $2;
        if (! exists $data_ref->{'gene'}->{$gene}->{'equiv'} ) {
            die "Cannot map gene $gene to an equivalence group\n";
        }
        my $equiv = $data_ref->{'gene'}->{$gene}->{'equiv'};
        print "$equiv\t$data\n";
    }
    else {
        die "From annotation file $annots, cannot parse: $annot\n";
    }
}
close $ANNOTS;
