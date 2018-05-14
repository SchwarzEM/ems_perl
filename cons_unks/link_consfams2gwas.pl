#!/usr/bin/env perl

use strict;
use warnings;
use List::MoreUtils qw(uniq);

my $family = $ARGV[0];
my $gwas   = $ARGV[1];

my $data_ref;

open my $GWAS, '<', $gwas or die "Can't open NIH GWAS index file $gwas: $!";
while (my $input = <$GWAS>) { 
    chomp $input;

# cat GWAS_OMIM/genome.gov.gwascatalog_18sep2013.txt | cut -f 8,15,25-26 | head --lines=5 ;
# 
# Disease/Trait   Mapped_gene     Context Intergenic
# Weight loss (gastric bypass surgery)    C5orf54 missense        0
# Weight loss (gastric bypass surgery)    MPPE1P - CAP2P1 Intergenic      1
# Weight loss (gastric bypass surgery)    AQP11   nearGene-5      0
# Weight loss (gastric bypass surgery)    SALL1 - UNGP1   Intergenic      1
# [...]
# Heart rate      RFX4;LOC100287944       intron;intron   0

    if ( $input =~ /\A (?: [^\t]* \t){7} ([^\t]*) \t (?: [^\t]* \t){6} ([^\t]*) \t (?: [^\t]* \t){9} [^\t]* \t ([^\t]*) /xms ) { 
        my $trait        = $1;
        my $mapped_genes = $2;
        my $intergenic   = $3;
        my @genes        = ();

        if ( ( $intergenic eq 1 ) and ( $mapped_genes =~ /\A (\S+) \s+ \- \s+ (\S+) \z/xms ) ) {
            my $gene1 = $1;
            my $gene2 = $2;
            @genes = ($gene1, $gene2);
            $trait = "intergenic_GWAS|\"$trait\"";
        }
        elsif ( ( $intergenic eq 0 ) and ( $mapped_genes =~ /\A (\S+) \z/xms ) ) { 
            my $gene_text = $1;
            @genes = split /;/, $gene_text;
            $trait = "intragenic_GWAS|\"$trait\"";
        }

        foreach my $hgnc_gene (@genes) { 
            if ( $hgnc_gene =~ /\A \S+ \z/xms ) { 
                $data_ref->{'hgnc_gene'}->{$hgnc_gene}->{'disease'}->{$trait} = 1;
            }
        }
    }
    else { 
        warn "From NIH GWAS index file $gwas, can't parse input: $input\n";
    }
}
close $GWAS or die "Can't close filehandle to NIH GWAS index file $gwas: $!";

open my $FAM, '<', $family or die "Can't open family table $family: $!";
while (my $input = <$FAM>) {
    chomp $input;
    my $hgnc_gene_gwas_text = q{};
    my @hgnc_gene_diseases  = ();

    # Sample input:
    # PF15189|DUF4582 13.0013 3 genes c_elegans|WBGene00012648|Y39A1A.9; human|ENSG00000180336|C17orf104; mouse|MGI:2686410|Gm1564 

    while ( $input =~ / human \| ENSG\d+ \| ([^;]+) ; /xmsg ) {
        my $hgnc_gene = $1;
        my @more_hgnc_gene_diseases = ();
        @more_hgnc_gene_diseases = map { "$hgnc_gene|$_" } sort keys %{ $data_ref->{'hgnc_gene'}->{$hgnc_gene}->{'disease'} }; 
        push @hgnc_gene_diseases, @more_hgnc_gene_diseases;
    }

    @hgnc_gene_diseases = sort @hgnc_gene_diseases ;
    @hgnc_gene_diseases = uniq @hgnc_gene_diseases ;

    if (@hgnc_gene_diseases) {
        $hgnc_gene_gwas_text = join '; ', @hgnc_gene_diseases;
    }

    print "$input\t$hgnc_gene_gwas_text\n";
}
close $FAM or die "Can't close filehandle to family table $family: $!";

