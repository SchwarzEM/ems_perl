#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

# Typical inputs:
# Entry	Gene names	Cross-reference (GeneTree)	Cross-reference (HGNC)
# A0A2U3TZM7	CACNB2	ENSGT00390000002740;	1402;
# A0A1W2PNQ7	ANKRD36	ENSGT00560000076605;	24079;
# F5H199	MSRB3	ENSGT00510000046802;	27375;
# G3V3U4	PSMA6	ENSGT00550000074807;	9535;
# 68431	HIST1H3A H3FA; HIST1H3B H3FL; HIST1H3C H3FC; HIST1H3D H3FB; HIST1H3E H3FD; HIST1H3F H3FI; HIST1H3G H3FH; HIST1H3H H3FK; HIST1H3I H3FF; HIST1H3J H3FJENSGT00760000118967;	4766;4776;4768;4767;4769;4773;4772;4775;4771;4774;
# Q8NG35	DEFB105A BD5 DEFB105 DEFB5; DEFB105B	ENSGT00390000002317;	18087;29930;
# Q5VWM5	PRAMEF9; PRAMEF15	ENSGT00760000119028;	26764;27996;
# P04745	AMY1A AMY1; AMY1B AMY1; AMY1C AMY1	ENSGT00390000002882;	474;475;476;
# P0C0S8	HIST1H2AG H2AFP; HIST1H2AI H2AFC; HIST1H2AK H2AFD; HIST1H2AL H2AFI; HIST1H2AM H2AFN	ENSGT00910000143981;	4737;4725;4726;4730;4735;

LOOP: while (my $input = <>) {
    chomp $input;
    my $uniprot_id    = q{};
    my $raw_gene_text = q{};
    my $raw_hgnc_text = q{};

    if ( $input =~ /\A (\S+) \t ([^\t]*) \t [^\t]* \t ( (?:\d+[;])+ ) \z/xms ) {
        $uniprot_id    = $1;
        $raw_gene_text = $2;
        $raw_hgnc_text = $3;
    }
    elsif ( $input =~ /\A (\S+) \t ([^\t]*) \t ( (?:ENSGT\d+[;])+ ) \t [^\t]* \z/xms ) {
        $uniprot_id    = $1;
        $raw_gene_text = $2;
        $raw_hgnc_text = $3;
    }
    else {
        warn "Cannot parse at all: $input\n";
        next LOOP;
    }

    if ( $uniprot_id and $raw_hgnc_text ) {
        my @orig_gene_texts = ();
        my @gene_texts      = ();
        my $gene_text_count = 0;

        @orig_gene_texts = split /[;]\s+/, $raw_gene_text ;
        foreach my $orig_gene_text (@orig_gene_texts) {
            if ( $orig_gene_text =~ /\A \s* (\S+) /xms ) {
                my $gene_text = $1;
                push @gene_texts, $gene_text;
            }
        }
        $gene_text_count = @gene_texts;

        my @hgnc_texts = ();
        my $hgnc_text_count = 0;
        @hgnc_texts = split /;/, $raw_hgnc_text;
        $hgnc_text_count = @gene_texts;

        if ( $gene_text_count != $hgnc_text_count ) {
            warn "Unequal counts of gene texts ($gene_text_count) and HGNC texts ($hgnc_text_count) in: $input\n";
        }

        my @full_gene_names = ();
        $hgnc_text_count--;

        foreach my $i (0..$hgnc_text_count) {
            my $hgnc_text = q{};
            my $gene_text = q{};
            $hgnc_text = $hgnc_texts[$i] if $hgnc_texts[$i]; 
            $gene_text = $gene_texts[$i] if $hgnc_texts[$i];
            if ( $hgnc_text =~ /\A \d+ \z/xms ) {
                $hgnc_text = "hgnc$hgnc_text";
            }
            if ($hgnc_text) {
                my $full_gene_name = "human|$hgnc_text";
                if ( $gene_text =~ /\A \S+ \z/xms ) { 
                    $full_gene_name = "human|$hgnc_text|$gene_text";
                }
                push @full_gene_names, $full_gene_name;
            }
        }

        my $full_gene_text = join '; ', @full_gene_names;
        $full_gene_text =~ s/[;]{2,}/;/g;
        print "$uniprot_id\t$full_gene_text\n";
    }
    else {
        warn "Cannot parse UniProt and gene IDs: $input\n";
    }
}
