#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $gene_type          = q{};
my $prelift            = q{};
my $postlift_badtrans  = q{};
my $postlift_no_ident  = q{};
my $postlift_all_ident = q{};

my $header = "Gene\tGene_type\tLiftover_status\n";

my @annot_files = ();

my $data_ref;

$gene_type = $ARGV[0] if $ARGV[0];

$prelift = $ARGV[1] if $ARGV[1];
push @annot_files, $prelift;

$postlift_badtrans = $ARGV[2] if $ARGV[2];
push @annot_files, $postlift_badtrans;

$postlift_no_ident = $ARGV[3] if $ARGV[3];
push @annot_files, $postlift_no_ident;

if ($ARGV[4]) {
    $postlift_all_ident = $ARGV[4];
    push @annot_files, $postlift_all_ident;
}
else {
    die "Format: tabulate_n2_lift_ncRNA.or.pseudo_annots_2018.08.06.01.pl\n",
        "    [gene-type annotation, e.g., \"ncRNA\" or \"pseudogene\"]\n",
        "    [at_least_prelift_gene_list]\n",
       	"    [at_least_postlift_badcoding_gene_list]\n",
       	"    [at_least_postlift_no_ident_gene_list]\n",
       	"    [at_least_postlift_all_ident_gene_list]\n",
        " > [annotation table]\n",
        ;
}

my @annot_texts = qw(
    no.lift
    lifted_bad.coding
    lifted_no.ident
    lifted_all.ident
);

my $INFILE;

foreach my $i (0..3) {
    my $annot_file = $annot_files[$i];
    my $annot_text = $annot_texts[$i];

    open $INFILE, '<', $annot_file;
    while (my $gene = <$INFILE>) {
        chomp $gene;
        $data_ref->{'gene'}->{$gene} = $annot_text;
    }
    close $INFILE;
}

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene (@genes) {
    print $header if $header;
    $header = q{};
    my $annot = $data_ref->{'gene'}->{$gene};
    print "$gene\t$gene_type\t$annot\n";
}

