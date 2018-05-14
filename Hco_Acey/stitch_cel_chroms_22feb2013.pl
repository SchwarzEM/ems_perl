#!/usr/bin/env perl

use strict;
use warnings;

my $chr_table = $ARGV[0];
my $hom_table = $ARGV[1];

my %gene2chr   = ( 'Strict_elegans_orthologs' => 'Strict_elegans_orth_chroms', );
my %gene2locus = ( 'Strict_elegans_orthologs' => 'Strict_elegans_orth_loci', );

open my $CHR, '<', $chr_table or die "Can't open chr. table $chr_table: $!";
while (my $input = <$CHR>) { 
    chomp $input;
    if ( $input =~ /\A (WBGene\d+) \t \S+ \t (\S+) \t (\S+) \t (\S+) \t (\S+) /xms ) {
        my $wbgene = $1;
        my $chr    = $2;
        my $ori    = $3;
        my $start  = $4;
        my $stop   = $5;
        my $annot  = $chr . q{:} . $start . q{-} . $stop . " [$ori]";
        $gene2chr{$wbgene}   = $chr;
        $gene2locus{$wbgene} = $annot;
    }
}
close $CHR or die "Can't close filehandle to chr. table $chr_table: $!";

open my $HOM, '<', $hom_table or die "Can't open homolog table $hom_table: $!";
while (my $input = <$HOM>) {
    chomp $input;
    if ( $input =~ /\A ([^\t]*) \t ([^\t]*) \t ([^\t]*) /xms ) { 
        my $gene         = $1;
        my $strict_orth  = $2;
        my $general_orth = $3;
        my $wb_gene_id   = $strict_orth;
        $wb_gene_id   =~ s/\A(WBGene\d+)\S*\z/$1/;
        print "$gene\t$strict_orth\t";
        if ( exists $gene2chr{$wb_gene_id} ) { 
            print "$gene2chr{$wb_gene_id}";
        }
        print "\t";
        if ( exists $gene2locus{$wb_gene_id} ) {
            print "$gene2locus{$wb_gene_id}";
        }
        print "\t$general_orth\n";
    }
}
close $HOM or die "Can't close filehandle to homolog table $hom_table: $!";

