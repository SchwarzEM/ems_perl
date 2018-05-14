#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::Util qw(max);
use List::MoreUtils qw(uniq);

my @genes = ();

my $data_ref;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) { 
        my $gene  = $1;
        my $value = $2;

        push @genes, $gene;
        $data_ref->{'gene'}->{$gene}->{'value'}->{$value} = 1;
    }
}

@genes = uniq(@genes);

foreach my $gene (@genes) {
    my @values    = sort { $a <=> $b } keys %{ $data_ref->{'gene'}->{$gene}->{'value'} };
    my $max_value = max(@values);
    print "$gene\t$max_value\n";
}

