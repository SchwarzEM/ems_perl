#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $gff  = q{};
my $list = q{};

$gff  = $ARGV[0] if $ARGV[0];
$list = $ARGV[1] if $ARGV[1];

my %banned = ();

if ( (! $gff ) or (! $list ) ) {
    die "Format: censor_agat_gff3_genes.pl [GFF3] [banned gene list] > [gene-line-censored GFF3]\n";
}

open my $LIST, '<', $list;
while ( my $input = <$LIST> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \z/xms ) {
        my $gene = $1;
        $banned{$gene} = 1;
    }
    else {
        die "From gene list $list, cannot parse: $input\n";
    }
}
close $LIST;

open my $GFF, '<', $gff;

while ( my $input = <$GFF> ) {
    chomp $input;
    if ( $input =~ /\A [^\t]* \t [^\t]* \t gene \t .+ \t ID=(\S+) /xms ) {
        my $gene = $1;
        if (! $banned{$gene} ) {
            print "$input\n";
        }
    }
    else {
        print "$input\n";
    }
}
close $GFF;
