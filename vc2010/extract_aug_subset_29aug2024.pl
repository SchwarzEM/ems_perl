#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $genelist = q{};
my $aug_gff  = q{};

$genelist = $ARGV[0] if $ARGV[0];
$aug_gff  = $ARGV[1] if $ARGV[1];

my %listed = ();

if ( (! $genelist ) or (! $aug_gff ) ) {
    die "Format: extract_aug_subset_29aug2024.pl [gene list] [AUGUSTUS GFF] > [gene list subset AUGUSTUS GFF]\n";
}

open my $GENELIST, '<', $genelist;
while (my $gene = <$GENELIST>) {
    chomp $gene;
    if ( $gene =~ /\A \S+ \z/xms ) {
        $listed{$gene} = 1;
    }
    else {
        die "From gene list $genelist, cannot parse: $gene\n";
    }
}
close $GENELIST;

open my $AUG_GFF, '<', $aug_gff;
while (my $input = <$AUG_GFF>) {
    chomp $input;
    if ( $input =~ /\A \S+ \t [^\t]* \t gene \t .* \t ID=(\S+) [^\t]* \z/xms ) {
        my $gene = $1;
        if ( exists $listed{$gene} ) {
            print "$input\n";
        }
    }
    elsif ( $input =~ /\A \S+ \t .* \t [^\t]* Parent= (\S+) \.t\d+ [^\t]* \z/xms ) {
        my $gene = $1;
        if ( exists $listed{$gene} ) {
            print "$input\n";
        }
    }
}
close $AUG_GFF;

