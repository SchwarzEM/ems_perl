#!/usr/bin/env perl

use strict;
use warnings;
use List::MoreUtils qw(uniq);

my $family        = $ARGV[0];
my $rnai_phenos   = $ARGV[1];
my $allele_phenos = $ARGV[2];

my $data_ref;

open my $RNAI, '<', $rnai_phenos or die "Can't open RNAi phenotype/no-phenotype data $rnai_phenos: $!";
while (my $input = <$RNAI>) { 
    chomp $input;

    # Sample input lines from cut -f 2,6,9 wb_data/ws235_gene2rnai_pheno_nopheno_18sep2013.tsv.txt ;
    # "WBGene00015175"	"embryonic lethal"
    # "WBGene00009334"                "embryonic lethal"

    # Get rid of xace's putting everything in "":
    $input =~ s/["]//g;

    if ( $input =~ /\A [^\t]* \t (WBGene\d+) \t (?: [^\t]* \t){3} ([^\t]*) \t (?: [^\t]* \t){2} ([^\t]*) /xms ) { 
        my $wbgene    = $1;
        my $phenotype = $2;
        my $no_pheno  = $3;

        if ( $phenotype =~ /\S/xms ) { 
            $data_ref->{'wb_gene'}->{$wbgene}->{'phenotype'}->{$phenotype} = 1;
        }
        if ( $no_pheno =~ /\S/xms ) {
            $data_ref->{'wb_gene'}->{$wbgene}->{'no_pheno'}->{$no_pheno} = 1;
        }
    }
    else { 
        warn "From RNAi phenotype/no-phenotype data $rnai_phenos, can't parse: $input\n";
    }
}
close $RNAI or die "Can't close filehandle to RNAi phenotype/no-phenotype data $rnai_phenos: $!";

open my $ALLELES, '<', $allele_phenos or die "Can't open allele phenotype data $allele_phenos: $!";
while (my $input = <$ALLELES>) {
    chomp $input;

    # Sample input from: cut -f 2,6 wb_data/ws235_gene2rnai_pheno_nopheno_18sep2013.tsv.txt | head --lines=3 
    # "WBGene00015175"	"embryonic lethal"
    # "WBGene00015175"	"embryonic lethal late emb"
    # "WBGene00015235"	"one cell arrest early emb"

    # Again, get rid of xace's putting everything in "":
    $input =~ s/["]//g;

    if ( $input =~ /\A [^\t]* \t (WBGene\d+) \t (?: [^\t]* \t){3} ([^\t]*) \t /xms ) { 
        my $wbgene    = $1;
        my $phenotype = $2;
        if ( $phenotype =~ /\S/xms ) {
            $data_ref->{'wb_gene'}->{$wbgene}->{'phenotype'}->{$phenotype} = 1;
        }
    }
    else { 
        warn "From allele phenotype data $allele_phenos, can't parse: $input\n";
    }
}
close $ALLELES or die "Can't close filehandle to allele phenotype data $allele_phenos: $!";

open my $FAM, '<', $family or die "Can't open family table $family: $!";
while (my $input = <$FAM>) {
    chomp $input;

    my $wb_gene_pheno_text    = q{};
    my $wb_gene_no_pheno_text = q{};

    my @wb_gene_phenos = ();
    my @wb_gene_no_phenos = ();

    while ( $input =~ /(WBGene\d+)/xmsg ) { 
        my $wbgene = $1;

        my @more_wb_gene_phenos = ();
        my @more_wb_gene_no_phenos = ();

        @more_wb_gene_phenos    = sort keys %{ $data_ref->{'wb_gene'}->{$wbgene}->{'phenotype'} };
        @more_wb_gene_no_phenos = sort keys %{ $data_ref->{'wb_gene'}->{$wbgene}->{'no_pheno'} };

        push @wb_gene_phenos, @more_wb_gene_phenos;
        push @wb_gene_no_phenos, @more_wb_gene_no_phenos;
    }

    @wb_gene_phenos = sort @wb_gene_phenos;
    @wb_gene_phenos = grep { /\S/ } uniq @wb_gene_phenos;

    @wb_gene_no_phenos = sort @wb_gene_no_phenos;
    @wb_gene_no_phenos = grep { /\S/ } uniq @wb_gene_no_phenos;

    if (@wb_gene_phenos) {
        $wb_gene_pheno_text = join '; ', @wb_gene_phenos;
        $wb_gene_pheno_text = "WBPheno: $wb_gene_pheno_text";
    }

    if (@wb_gene_no_phenos) {
        $wb_gene_no_pheno_text = join '; ', @wb_gene_no_phenos;
        $wb_gene_no_pheno_text =  "WBNonPheno: $wb_gene_no_pheno_text";
    }       
    print "$input\t$wb_gene_pheno_text\t$wb_gene_no_pheno_text\n";
}
close $FAM or die "Can't close filehandle to family table $family: $!";

