#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    # require at least one TAIR ID in the gene info
    if ( $input =~ /\A (\S+) \t ([^\t]* At[0-9A-Za-z]g\d+ [^\t]*) \z/xms ) {
        my $uniprot_id    = $1;
        my $raw_gene_text = $2;

        # Note: splitting on ';' alone does not work, because Arabidopsis has idiotic gene names like 'PHT3;3'
        # Also note: UniProt has junk entries like "At1gXXXXX F13E17.2", but the grep should filter them out.
        my @gene_infos    = grep { $_ =~ /At[0-9A-Za-z]g\d+/ } split /[;][ ]/, $raw_gene_text;
        my @full_genes    = ();

        foreach my $gene_info (@gene_infos) {
            $gene_info =~ s/\A\s+//;
            if ( $gene_info =~ /\A ([^\t]*?) (At[0-9A-Za-z]g\d+)/xms ) { 
                my $gene_name = $1;
                my $tair_id   = $2;
                my $full_name = "arabidopsis|$tair_id";
                if ( $gene_name =~ /\A (\S+)/xms ) {
                     $gene_name = $1;
                     $full_name = "arabidopsis|$tair_id|$gene_name";
                }
                push @full_genes, $full_name;
            }
            else { 
                die "Cannot parse gene info ($gene_info) in: $input\n";
            }
        }
        @full_genes = sort @full_genes;
        my $full_gene_text = join '; ', @full_genes;
        print "$uniprot_id\t$full_gene_text\n";
    }
    else {
        warn "Cannot parse at all: $input\n";
    }
}

