#!/usr/bin/env perl

use strict;
use warnings;
use List::MoreUtils qw(uniq);

my $family    = $ARGV[0];
my $phenos    = $ARGV[1];
my $no_phenos = $ARGV[2];

my $data_ref;

open my $PHENOS, '<', $phenos or die "Can't open phenotype data $phenos: $!";
while (my $input = <$PHENOS>) { 
    chomp $input;

    # Sample input lines from cut -f 1,5-6 wb_data/wormmine_phenos_17sep2013.tsv.txt | head ;
    # Gene > WormBase Gene ID	Gene > Alleles > Phenotypes Observed > Identifier	Gene > Alleles > Phenotypes Observed > Name
    # WBGene00000001	WBPhenotype:0000964	DMPP resistant
    # WBGene00000001	WBPhenotype:0000674	slow development

    if ( $input =~ /\A (WBGene\d+) \t (?: [^\t]* \t){3} WBPhenotype:\d+ \t ([^\t]+) \z/xms ) { 
        my $wbgene    = $1;
        my $phenotype = $2;
        $data_ref->{'wb_gene'}->{$wbgene}->{'phenotype'}->{$phenotype} = 1;
    }
}
close $PHENOS or die "Can't close filehandle to phenotype data $phenos: $!";

open my $NO_PHENOS, '<', $no_phenos or die "Can't open no-phenotype data $no_phenos: $!";
while (my $input = <$NO_PHENOS>) {
    chomp $input;
    # Input as before!
    if ( $input =~ /\A (WBGene\d+) \t (?: [^\t]* \t){3} WBPhenotype:\d+ \t ([^\t]+) \z/xms ) {
        my $wbgene    = $1;
        my $no_phenotype = $2;
        $data_ref->{'wb_gene'}->{$wbgene}->{'no_phenotype'}->{$no_phenotype} = 1;
    }
}
close $NO_PHENOS or die "Can't close filehandle to no-phenotype data $no_phenos: $!";

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
        @more_wb_gene_no_phenos = sort keys %{ $data_ref->{'wb_gene'}->{$wbgene}->{'no_phenotype'} };

        push @wb_gene_phenos, @more_wb_gene_phenos;
        push @wb_gene_no_phenos, @more_wb_gene_no_phenos;
    }

    @wb_gene_phenos = sort @wb_gene_phenos;
    @wb_gene_phenos = uniq @wb_gene_phenos;

    @wb_gene_no_phenos = sort @wb_gene_no_phenos;
    @wb_gene_no_phenos = uniq @wb_gene_no_phenos;

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
