#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\S/xms ) {
        $data_ref->{'annot'}->{$input}->{'count'}++
    }
}

my @annots = sort keys %{ $data_ref->{'annot'} };
foreach my $annot (@annots) {
    my $count = $data_ref->{'annot'}->{$annot}->{'count'};
    print "$annot\t$count\n";
}




