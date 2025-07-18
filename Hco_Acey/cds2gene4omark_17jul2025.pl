#!/usr/bin/env perl

use strict;
use autodie;
use warnings;

my $data_ref;

while ( my $input = <> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) /xms ) {
        my $cds  = $1;
        my $gene = $2;
        $data_ref->{'gene'}->{$gene}->{'cds'}->{$cds} = 1;
    }
    else {
       die "Cannot parse input: $input\n";
    }
}

my @genes = sort keys %{ $data_ref->{'gene'} };

foreach my $gene (@genes) {
    my @cdses = sort keys %{ $data_ref->{'gene'}->{$gene}->{'cds'} };
    my $cds_list = join ';', @cdses;
    print "$cds_list\n";
}

