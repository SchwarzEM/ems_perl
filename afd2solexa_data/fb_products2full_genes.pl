#!/usr/bin/env perl

# ens_products2full_genes.pl -- Erich Schwarz <ems394@cornell.edu>, 3/21/2020.
# Purpose: given Ensembl proteome headers for human or mouse, build a ENSG \t pubname table which ens2fullnames.pl can then use.

use strict;
use warnings;

my $gene_id     = q{};
my $gene_symbol = q{};

my $data_ref;

# sample input line:
# name=nAChRalpha2-PA; parent=FBgn0000039,[...]
# name=a-PA; parent=FBgn0000008,FBtr0071763;

while (my $input = <>) {
    if ( $input =~ /\A > /xms ) { 
        if ( $input !~ / \A > .+ \s name=\S+ [-]P[A-Z]+ [;] \s parent=FBgn\d+[,] /xms ) { 
            die "Can't parse FASTA header line: $input\n";
        }
        if ( $input =~ / \A > .+ \s name=(\S+) [-]P[A-Z]+ [;] \s parent=(FBgn\d+)[,] /xms ) {
            $gene_symbol = $1;
            $gene_id     = $2;
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

