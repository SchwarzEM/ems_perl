#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $uprot2gene = q{};
my $gene_annot = q{};

my $data_ref;

my @annotated_genes = ();

$uprot2gene = $ARGV[0] if $ARGV[0];
$gene_annot = $ARGV[1] if $ARGV[1];

if ( (! $uprot2gene ) or (! $gene_annot ) ) {
     die "Format: add_uniprots2gene_annots_30sep2024.pl [UniProt/GenBank/Gene] [Gene annots] > [Gene/UniProt/GenBank/annots] ;\n";
}

open my $GENE_ANNOT, '<', $gene_annot;
while (my $input = <$GENE_ANNOT>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (.*) \z/xms ) {
        # Unlike link_uprots2genes_05may2024.pl, we record only the actual annotation text and do *not* include a leading gene name.
        my $gene  = $1;
        my $annot = $2;
        # We actually accept 'Gene' as a gene ID as a hack to copy over all of our header columns!
        if ( exists $data_ref->{'gene'}->{$gene}->{'annot'} ) {
            die "Redundant annotations of gene $gene\n"
        }
        else {
            # I'm making a decision here to allow genes without any annotation text to be added, as long as they at least have \t:
            if ( $annot !~ /\S /xms ) {
                warn "For gene $gene in annotation table $gene_annot, empty annotation text: \"$annot\"\n";
            }
            if ( $annot !~ /\t /xms ) {
                die "For gene $gene in annotation table $gene_annot, annotation text has no tabs: \"$annot\"\n";
            }

            # If a gene makes it past those filters, add it to the annotation set:
            push @annotated_genes, $gene;
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

    # Sample input:
    # UniProt_ID      GenBank_ID      Gene
    # A0AA39GMM2      QR680_019442    Sherm_2022.07.24.54.g339
    # A0AA39GMN1      QR680_019434    Sherm_2022.07.24.39.g330

    # We are mapping *genes* to *protein IDs*, so we must accomodate cases where 1 gene encodes 2+ proteins.
    if ( $input =~ /\A (\S+) \t (\S+) \t (\S+) \z/xms ) {
        my $uprot_id = $1;
        my $gbank_id = $2;
        my $gene     = $3;
        if (! exists $data_ref->{'gene'}->{$gene}->{'annot'} ) {
            warn "No annotation of gene $gene\n"
        }
        else {
            my $uprot_etc = "$uprot_id|$gbank_id";
            $data_ref->{'gene'}->{$gene}->{'uprot_etc'}->{$uprot_etc} = 1;
        }
    }
    else {
        die "From gene annotation file $gene_annot, cannot parse: $input\n";
    }
}
close $UPROT2GENE;

foreach my $annotated_gene (@annotated_genes) {
    my $annot      = $data_ref->{'gene'}->{$annotated_gene}->{'annot'};
    my $uprot_text = "None|None";
    if ( exists $data_ref->{'gene'}->{$annotated_gene}->{'uprot_etc'} ) {
        my @uprot_etcs = sort keys %{ $data_ref->{'gene'}->{$annotated_gene}->{'uprot_etc'} };
        $uprot_text    = join '; ', @uprot_etcs;
    }
    print "$annotated_gene\t$uprot_text\t$annot\n";
}

