#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

while ( my $input = <> ) {
    chomp $input;
    if ( $input =~ /\A ( (\S+) \.t\d+) \z/xms ) {
        my $tx   = $1;
        my $gene = $2;
        $data_ref->{'gene'}->{$gene}->{'tx'}->{$tx} = 1;
    }
    else {
        die "Cannot parse: $input\n";
    }
}

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene (@genes) {
    my @txs = sort keys %{ $data_ref->{'gene'}->{$gene}->{'tx'} };
    my $tx_text = join '; ', @txs;
    print "$gene\t$tx_text\n";
}
