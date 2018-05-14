#!/usr/bin/env perl

use strict;
use warnings;

my $pheno_table = $ARGV[0];
my $drug_table  = $ARGV[1];

my %gene2pheno = ();

open my $PHENOS, '<', $pheno_table or die "Can't open pheno table $pheno_table: $!";
while (my $input = <$PHENOS>) {
    chomp $input;
    if ( $input =~ /\A (WBGene\d+\S+) \t ([^\t]+) /xms ) { 
        my $gene  = $1;
        my $pheno = $2;
        $pheno =~ s/\"//g;
        $pheno =~ s/WBPhenotype: //;
        $pheno = ucfirst $pheno;
        $gene2pheno{$gene} = $pheno;
    }
}
close $PHENOS or die "Can't close filehandle to pheno table $pheno_table: $!";

open my $DRUGS, '<', $drug_table or die "Can't open drug target table $drug_table: $!";
while (my $input = <$DRUGS>) {
    chomp $input;
    if ( $input =~ /\A .+ \t (WBGene\d+\S+) \z/xms ) { 
        my $gene  = $1;
        my $pheno = q{};
        if ( exists $gene2pheno{$gene} ) { 
            $pheno = $gene2pheno{$gene};
        }
        print "$input\t$pheno\n";
    }
}
close $DRUGS or die "Can't close filehandle to drug target table $drug_table: $!";

