#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

# Typical inputs:
# A1Z8J4	nompA CG13207 Dmel_CG13207	FBgn0016047;
# Q7KN74	l(2)k14710 CG8325 Dmel_CG8325	FBgn0021847;
# Q7K7B0	Ugt37c1 ec EG:EG0003.4 Ugt37c1-RA CG8652 Dmel_CG8652	FBgn0026754;
# Q494G1	CG9864 Dmel_CG9864	FBgn0034490;
# A1ZB48	CG14500 Dmel_CG14500	FBgn0034318;
# Q9W2N8	CG10543 Dmel_CG10543	FBgn0034570;
# Q7K3W2	CG8728-RA CG8728 Dmel_CG8728	FBgn0033235;
# Q5BIB8	CG17443 CG41520 Dmel_CG41520	FBgn0087011;
# Q7K2J4	dmpd CG11866 Dmel_CG11866	FBgn0033486;

while (my $input = <>) {
    chomp $input;
    # require at least one TAIR ID in the gene info
    if ( $input =~ /\A (\S+) \t ([^\t]* CG\d+ [^\t]*) \t ([^\t]* FBgn\d+[;] [^\t]*)\z/xms ) {
        my $uniprot_id    = $1;
        my $raw_gene_text = $2;
        my $raw_fb_text   = $3;

        my @gene_texts = ();
        my $gene_text_count = 0;

        # Clean out all the junky derivatives of CG\d+:
        $raw_gene_text =~ s/\S+CG\d+//g;

        if ( $raw_gene_text =~ / \b (CG\d+) \b /xms ) {
            while ($raw_gene_text =~ /\A ([^\t]*?) \b (CG\d+) \b /xmsg ) { 
                my $gene_info  = $1;
                my $flybase_id = $2;

                if ( $gene_info =~ /\A \s* (\S+) /xms ) {
                    $flybase_id = $1;
                }
                push @gene_texts, $flybase_id;
            }
        }

        $gene_text_count = @gene_texts;
 
        my @fb_texts = ();
        my $fb_text_count = q{};

        $raw_fb_text =~ s/[;]\z//;
        @fb_texts = split /[;]/, $raw_fb_text;
        $fb_text_count = @fb_texts;

        if ( $gene_text_count != $fb_text_count ) {
            warn "Unequal counts of gene texts ($gene_text_count) and FlyBase texts ($fb_text_count) in: $input\n";
        }

        my @full_gene_names = ();
        $fb_text_count--;

        foreach my $i (0..$fb_text_count) {
            my $fb_text   = $fb_texts[$i]; 
            my $gene_text = $gene_texts[$i];
            my $full_gene_name = "drosophila|$fb_text";
            if ( $gene_text =~ /\A \S+ \z/xms ) { 
                $full_gene_name = "drosophila|$fb_text|$gene_text";
            }
            push @full_gene_names, $full_gene_name;
        }

        my $full_gene_text = join '; ', @full_gene_names;
        $full_gene_text =~ s/[;]{2,}/;/g;
        print "$uniprot_id\t$full_gene_text\n";
    }
    else {
        warn "Cannot parse at all: $input\n";
    }
}

