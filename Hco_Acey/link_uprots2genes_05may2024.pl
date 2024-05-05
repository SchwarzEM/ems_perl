#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $uprot2gene = q{};
my $gene_annot = q{};

my $data_ref;

$uprot2gene = $ARGV[0] if $ARGV[0];
$gene_annot = $ARGV[1] if $ARGV[1];

if ( (! $uprot2gene ) or (! $gene_annot ) ) {
     die "Format: link_uprots2genes_05may2024.pl [UniProt/GenBank/Gene] [Gene annots] > [UniProt/GenBank/Gene/annots] ;\n";
}

open my $GENE_ANNOT, '<', $gene_annot;
while (my $input = <$GENE_ANNOT>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t .* \z/xms ) {
        my $gene = $1;
        # We actually accept 'Gene' as a gene ID as a hack to copy over all of our header columns!
        if ( exists $data_ref->{'gene'}->{$gene}->{'annot'} ) {
            die "Redundant annotations of gene $gene\n"
        }
        else {
            $data_ref->{'gene'}->{$gene}->{'annot'} = $input;
        }
    }
    else {
        die "From gene annotation file $gene_annot, cannot parse: $input\n";
    }
}
close $GENE_ANNOT;

open my $UPROT2GENE, '<', $uprot2gene;
while (my $input = <$UPROT2GENE>) {
    chomp $input;
    if ( $input =~ /\A (\S+ \t \S+) \t (\S+) \z/xms ) {
        my $uprot_etc = $1;
        my $gene      = $2;
        if (! exists $data_ref->{'gene'}->{$gene}->{'annot'} ) {
            die "No annotation of gene $gene\n"
        }
        else {
            my $annot =  $data_ref->{'gene'}->{$gene}->{'annot'};
            print "$uprot_etc\t$annot\n";
        }
    }
    else {
        die "From gene annotation file $gene_annot, cannot parse: $input\n";
    }
}
close $UPROT2GENE;
