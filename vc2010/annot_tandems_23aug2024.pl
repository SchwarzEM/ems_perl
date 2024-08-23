#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my $header = "Gene\tTandem_copies";

while ( my $gene = <> ) {
    chomp $gene;
    if ( $gene =~ /\A (\S+)_\d+ \z/xms ) {
        my $parent_gene = $1;
        $data_ref->{'parent_gene'}->{$parent_gene}->{'gene'}->{$gene} = 1;
    }
}

my @parents = sort keys %{ $data_ref->{'parent_gene'} };
foreach my $parent (@parents) {
    my @tandems = sort keys %{ $data_ref->{'parent_gene'}->{$parent}->{'gene'} };
    my $count = @tandems;
    $count++;

    print "$header\n" if $header;
    $header = q{};

    print "$parent\tTandem_copies[x$count]\n";
}
