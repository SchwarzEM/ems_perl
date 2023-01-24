#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $tx2gene_tsv = q{};
my $tx_list = q{};

my %tx2gene    = ();
my %seen_txs   = ();
my %seen_genes = ();

$tx2gene_tsv = $ARGV[0] if $ARGV[0];
$tx_list     = $ARGV[1] if $ARGV[1];

if ( (! -r $tx2gene_tsv ) or (! -r $tx_list ) ) {
    die "Format: map_tx_to_genes.pl [transcript to gene TSV] [transcript list] > [sorted unique gene list]\n";
}

open my $TX2GENE, '<', $tx2gene_tsv;
while ( my $input = <$TX2GENE> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $tx   = $1;
        my $gene = $2;
        if ( exists $seen_txs{$tx} ) {
            die "Cannot map transcript $tx to a gene 2+ times\n";
        }
        $tx2gene{$tx}  = $gene;
        $seen_txs{$tx} = 1;
    }
    else {
        die "From transcript to gene TSV $tx2gene_tsv, cannot parse: $input\n";
    }
}
close $TX2GENE;

open my $TX_LIST, '<', $tx_list;
while ( my $input = <$TX_LIST> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \z/xms ) {
        my $tx   = $1;
        if (! exists $tx2gene{$tx} ) {
            warn "Cannot map transcript $tx to any gene\n";
        }
        else {
            my $gene = $tx2gene{$tx};
            $seen_genes{$gene} = 1;
        }
    }
    else {
        die "From transcript list $tx_list, cannot parse: $input\n";
    }
}
close $TX_LIST;

my @gene_list = sort keys %seen_genes;
foreach my $gene (@gene_list) {
    print "$gene\n";
}
