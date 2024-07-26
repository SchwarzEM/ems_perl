#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

while (my $input = <>) {
    chomp $input;

    if ( $input =~ /\A \S+ \t [^\t]* \t gene \t .* \t gene_id \s+ \"(WBGene\d+) \" .* gene_name \s+ \"(\S+)\" /xms ) {
        my $gene_id   = $1;
        my $gene_name = $2;
        $data_ref->{'gene_id'}->{$gene_id}->{'gene_name'} = $gene_name;
    }
    elsif ( $input =~ /\A \S+ \t [^\t]* \t gene \t /xms ) {
        die "Cannot parse gene id/name input: $input\n";
    }

    if ( $input =~ /\A \S+ \t [^\t]* \t transcript \t .* \t gene_id \s+ \"(WBGene\d+) \" .* transcript_id \s+ \"(\S+)\" /xms ) {
        my $gene_id   = $1;
        my $transcript_id = $2;
        $data_ref->{'gene_id'}->{$gene_id}->{'transcript_id'} = $transcript_id;
    }
    elsif ( $input =~ /\A \S+ \t [^\t]* \t transcript \t /xms ) {
        die "Cannot parse gene id/name input: $input\n";
    }
}

my @gene_ids = sort keys %{ $data_ref->{'gene_id'} };

foreach my $gene_id (@gene_ids) {
    my $gene_name     = $data_ref->{'gene_id'}->{$gene_id}->{'gene_name'};
    my $transcript_id = $data_ref->{'gene_id'}->{$gene_id}->{'transcript_id'};
    print "$gene_id\t$transcript_id\t$gene_name\n";
}
