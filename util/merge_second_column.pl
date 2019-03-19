#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $index = $1;
        my $annot = $2;
        $data_ref->{'index'}->{$index}->{'annot'}->{$annot} = 1;
    }
    else { 
        die "Cannot parse input line: $input\n";
    }
}

my @indexes = sort keys %{ $data_ref->{'index'} };

foreach my $index (@indexes) {
    my @annots     = sort keys %{ $data_ref->{'index'}->{$index}->{'annot'} };
    my $annot_text = join '; ', @annots;
    print "$index\t$annot_text\n";
}


