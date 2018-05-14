#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (TF\d+) \t (\S+) \t (\S+) \z/xms ) { 
        my $treefam = $1;
        my $species = $2;
        my $gene    = $3;
        if ( $species eq 'homo_sapiens' ) { 
            $data_ref->{'hum_treefam'}->{$treefam} = 1;
        }
        elsif ( $species eq 'caenorhabditis_elegans' ) {
            $data_ref->{'worm_treefam'}->{$treefam}->{'gene'}->{$gene} = 1;
        }
    }
    else { 
        die "Can't parse input: $input\n";
    }
}

my @ok_treefams = grep { exists $data_ref->{'worm_treefam'}->{$_} } sort keys %{ $data_ref->{'hum_treefam'} };

foreach my $ok_treefam (@ok_treefams) { 
    my @genes = sort keys %{ $data_ref->{'worm_treefam'}->{$ok_treefam}->{'gene'} };
    foreach my $gene (@genes) {
        print "$gene\t$ok_treefam\n";
    }
}

