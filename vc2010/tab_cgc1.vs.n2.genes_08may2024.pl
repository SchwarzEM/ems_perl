#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = "Gene\tMapped_N2_gene";

my $data_ref;

while (my $input = <>) {
    chomp $input;
    # ID=chrI.g1.t1.cds;Parent=chrI.g1.t1     WBGene00022277
    if ( $input =~ /\A ID=\S+\.g\d+\.t\d+\.cds;Parent=(\S+\.g\d+)\.t\d+ \t (WBGene\d+) \z/xms ) {
        my $cgc1_gene = $1;
        my $n2_gene   = $2;
        $data_ref->{'cgc1_gene'}->{$cgc1_gene}->{'n2_gene'}->{$n2_gene} = 1;
    }
    elsif ( $input !~ /\A \S+ \t \. \z/xms ) {
        die "Cannot parse input: $input\n";
    }
}

my @cgc1_genes = sort keys %{ $data_ref->{'cgc1_gene'} };
foreach my $cgc1_gene (@cgc1_genes) {
    my @n2_genes = sort keys %{ $data_ref->{'cgc1_gene'}->{$cgc1_gene}->{'n2_gene'} };
    my $n2_text = join '; ', @n2_genes;

    print "$header\n" if $header;
    $header = q{};

    print "$cgc1_gene\t$n2_text\n";
}
