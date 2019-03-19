#!/usr/bin/env perl

# wormpep2uniprot_wbgene_table.pl -- Erich Schwarz <ems394@cornell.edu>, 11/2/2018.
# Purpose: given headers from wormpep (in format seen for WS267, 11/2018) emit "[uniprot] \t worm | WBGene | [human-readable]" TSV table.

use strict;
use warnings;
use autodie ;

my $data_ref;

# Typical wormpep WS267 lines:
# >2L52.1b wormpep=CE50569 gene=WBGene00007063 status=Confirmed uniprot=A0A0K3AWR5 insdc=CTQ86426.1
# >2RSSE.1a wormpep=CE32785 gene=WBGene00007064 locus=rga-9 status=Partially_confirmed uniprot=A4F337 insdc=CCD61138.1

while (my $input = <>) { 
    chomp $input;
    if ($input =~ /\A [>] (\S+) \b .+ (WBGene\d+) .+ uniprot=(\S+) /xms) { 
        my $gene        = $1;
        my $wbgene      = $2;
        my $uniprot     = $3;
        my $locus       = q{};

        $gene =~ s/[a-z]\z//;

        if ( $input =~ / locus = (\S+) /xms ) {
            $gene = $1;
        }

        my $full_gene = "elegans|$wbgene|$gene";
        $data_ref->{'uniprot'}->{$uniprot}->{'full_gene'}->{$full_gene} = 1;
    }
    elsif ( $input =~ /\A [>] /xms ) { 
        warn "Cannot parse FASTA header line: $input\n";
    }
}

my @uniprots = sort keys %{ $data_ref->{'uniprot'} };
foreach my $uniprot (@uniprots) {
    my @full_genes = sort keys %{ $data_ref->{'uniprot'}->{$uniprot}->{'full_gene'} };
    my $full_gene_text = join '; ', @full_genes;
    print "$uniprot\t$full_gene_text\n";
}

