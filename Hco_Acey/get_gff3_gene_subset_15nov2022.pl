#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $genes = q{};
my $gff3  = q{};

my %ok_genes = ();

$genes = $ARGV[0] if $ARGV[0];
$gff3  = $ARGV[1] if $ARGV[1];

if ( (! -e $genes ) or (! -e $gff3 ) ) {
    die "Format: get_gff3_gene_subset_15nov2022.pl [selected gene list] [GFF3 with gene-annotated lines] > [GFF3 subset with only selected genes\n";
}

open my $GENES, '<', $genes;
while (my $input = <$GENES>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \z/xms ) {
        my $gene = $1;
        $ok_genes{$gene} = 1;
    }
    else {
        die "In gene list $genes, cannot parse: $input\n";
    }
}
close $GENES;

open my $GFF3, '<', $gff3;
while (my $input = <$GFF3>) {
    chomp $input;
    if ( $input =~ / gene=([^\s;]+) /xms ) {
        my $gene = $1;
        if ( exists $ok_genes{$gene} ) {
            print "$input\n";
        }
    }
    else {
        die "In GFF3 $gff3, cannot parse: $input\n";
    }
}
close $GFF3;
