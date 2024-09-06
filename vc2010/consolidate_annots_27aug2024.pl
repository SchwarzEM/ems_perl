#!/usr/bin/env perl

use strict;
use autodie;
use strict;

use List::MoreUtils qw(uniq);

my @genes = ();
my $data_ref;

while (my $input = <> ) {
    chomp $input;
    if ( $input =~ /\A ([^\t]+) \t ([^\t]+) \z/xms ) {
        my $gene  = $1;
        my $annot = $2;
        push @genes, $gene;
        $data_ref->{'gene'}->{$gene}->{'annot'}->{$annot} = 1;
    }
    else {
        die "Cannot parse input: $input\n";
    }
}

@genes = uniq(@genes);

foreach my $gene (@genes) {
    if ( exists $data_ref->{'gene'}->{$gene}->{'annot'} ) {
        my @annots = sort keys %{ $data_ref->{'gene'}->{$gene}->{'annot'} };
        my $annot_text = join '; ', @annots;
        print "$gene\t$annot_text\n";
    }
}

