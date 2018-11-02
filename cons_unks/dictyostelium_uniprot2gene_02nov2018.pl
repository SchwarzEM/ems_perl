#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

# Sample input:
# Q556I0	DDB0202985 DDB0217436	DDB_G0272879;DDB_G0274047;
# Q556T9	DDB0203032 DDB0217327	DDB_G0273279;DDB_G0273827;
# Q556M2	DDB0203006 DDB0217390	DDB_G0272901;DDB_G0273955;
# O15916	rsc22 DDB0214884	DDB_G0279795;
# Q8MML4	rcdBB DDB0185118	DDB_G0274551;
# Q556P0	DDB0168013 DDB0217386	DDB_G0272602;DDB_G0273947;
# Q557C2	DDB0168135 DDB_G0273661	DDB_G0273241;DDB_G0273661;
# Q556N4	rsc11-2 rsc11 DDB0185119 DDB_G0295699	DDB_G0272791;DDB_G0295699;
# Q54W23	drpp40 DDB0205891	DDB_G0279961;
# Q54P86	TPS6 DDB0186155	DDB_G0284707;

while (my $input = <>) {
    chomp $input;
    # require at least one Dicty ID in the gene and dicytBase info
    if ( $input =~ /\A (\S+) \t ([^\t]* DDB(?:_){0,1}(?:G){0,1}\d+ [^\t]*) \t ([^\t]* DDB_G\d+[;] [^\t]*) \z/xms ) {
        my $uniprot_id     = $1;
        my $raw_gene_text  = $2;
        my $raw_dicty_text = $3;

        my @gene_texts = ();
        while ($raw_gene_text =~ / ([^\t]*?) (DDB(?:_){0,1}(?:G){0,1}\d+) /xmsg ) { 
            my $gene_info = $1;
            my $dicty_id  = $2;

            if ( $gene_info =~ /\A (\S+)/xms ) {
                $dicty_id = $1;
            }
            
            push @gene_texts, $dicty_id;
        }
        my $gene_text_count = @gene_texts;

        my @dicty_texts = ();
        my $dicty_text_count = q{};

        $raw_dicty_text =~ s/[;]\z//;
        @dicty_texts = split /[;]/, $raw_dicty_text;
        $dicty_text_count = @dicty_texts;

        if ( $gene_text_count != $dicty_text_count ) {
            warn "Unequal counts of gene texts ($gene_text_count) and dicty texts ($dicty_text_count) in: $input\n";
        }

        my @full_gene_names = ();
        $dicty_text_count--;
        my @i_series = (0..$dicty_text_count);
        @i_series = reverse @i_series;

        foreach my $i (@i_series) {
            my $dicty_text = $dicty_texts[$i]; 
            my $gene_text  = $gene_texts[$i];
            my $full_gene_name = "dicty|$dicty_text|$gene_text";
            push @full_gene_names, $full_gene_name;
        }

        my $full_gene_text = join '; ', @full_gene_names;
        print "$uniprot_id\t$full_gene_text\n";
    }
    else {
        warn "Cannot parse at all: $input\n";
    }
}

