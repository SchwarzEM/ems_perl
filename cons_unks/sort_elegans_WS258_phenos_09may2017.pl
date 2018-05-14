#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \t (\S+) \t (.+) \z/xms ) { 
        my $wbgene    = $1;
        my $cgc       = $2;
        my $cds       = $3;
        my $phenotext = $4;

        my $gene = q{};

        if ( $cgc eq $cds ) {
            $gene = "$wbgene|$cds";
        }
        else {
            $gene = "$wbgene|$cds|$cgc";
        }

        $phenotext =~ s/N\.A\.//g; 
        $phenotext =~ s/\t/,/g;

        $phenotext =~ s/[,]\z//;

        my @phenos = split /[,]/, $phenotext;
        @phenos    = sort @phenos;

        foreach my $pheno (@phenos) {
            $data_ref->{'gene'}->{$gene}->{'pheno'}->{$pheno} = 1;
        }
    }
    else {
       die "Cannot parse input: $input\n";
    }
}

my @genes = sort keys %{ $data_ref->{'gene'} };

foreach my $gene (@genes) {
    my @phenos = sort keys %{ $data_ref->{'gene'}->{$gene}->{'pheno'} };
    foreach my $pheno (@phenos) {
        print "$gene\t$pheno\n";
    }
}


