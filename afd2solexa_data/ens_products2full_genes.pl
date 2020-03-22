#!/usr/bin/env perl

# ens_products2full_genes.pl -- Erich Schwarz <ems394@cornell.edu>, 3/21/2020.
# Purpose: given Ensembl proteome headers for human or mouse, build a ENSG \t pubname table which ens2fullnames.pl can then use.

use strict;
use warnings;

my $gene_id     = q{};
my $gene_symbol = q{};

my $data_ref;

# sample input line:
# >ENSP00000362111.4 pep chromosome:GRCh38:X:100627108:100636806:-1 \
# gene:ENSG00000000003.15 transcript:ENST00000373020.9 gene_biotype:protein_coding transcript_biotype:protein_coding \
# gene_symbol:TSPAN6 description:tetraspanin 6 [Source:HGNC Symbol;Acc:HGNC:11858]

while (my $input = <>) {
    if ( $input =~ /\A > /xms ) { 
        if ( $input !~ / \A > .+ \s gene:\S+ \s .+ \s gene_symbol:\S+ \s /xms ) { 
            die "Can't parse FASTA header line: $input\n";
        }
        if ( $input =~ / \A > .+ \s gene:(\S+) \s .+ \s gene_symbol:(\S+) \s /xms ) {
            $gene_id     = $1;
            $gene_symbol = $2;
            # Note: ENSEMBL allows two genes to share a single gene symbol, unlike WormBase, so don't ban that.
            $data_ref->{'gene_id'}->{$gene_id}->{'gene_symbol'} = $gene_symbol;
       }
   }
}

my @gene_ids = sort keys %{ $data_ref->{'gene_id'} };

foreach my $gene_id1 (@gene_ids ) {
    my $gene_symbol1 = $data_ref->{'gene_id'}->{$gene_id1}->{'gene_symbol'};
    print "$gene_id1\t$gene_symbol1\n";
}

