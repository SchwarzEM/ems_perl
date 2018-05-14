#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $bad_genes  = $ARGV[0];
my $gene2annot = $ARGV[1];

my $data_ref;

if ( (! $bad_genes) or (! $gene2annot) ) { 
    die "Format: censor_gene2annot_22feb2014.pl [bad genes list] [previous gene-to-annotation list (e.g., gene-to-locus)] > [new, censored gene-to-annotation list]\n";
}

open my $BAD_GENES, '<', $bad_genes;
while (my $input = <$BAD_GENES>) {
    chomp $input;
    if ( $input =~ /\A \S+ \z/xms ) {
        $data_ref->{'bad_gene'}->{$input} = 1;
    }
    else { 
        die "From bad genes file $bad_genes, can't parse: $input\n";
    }
}
close $BAD_GENES;

open my $GENE2ANNOT, '<', $gene2annot;
while (my $input = <$GENE2ANNOT>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t \S+ .* \z/xms ) { 
        my $gene = $1;
        if (! exists $data_ref->{'bad_gene'}->{$gene} ) {
            print "$input\n";
        }
    }
    else { 
        die "From gene-to-annotation file $gene2annot, can't parse: $input\n";
    }
}
close $GENE2ANNOT;

