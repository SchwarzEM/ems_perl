#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) { 
        my $gene    = $1;
        my $treefam = $2;
        $data_ref->{'gene'}->{$gene}->{'treefam'}->{$treefam} = 1;
    }
    else { 
        die "Can't parse input: $input\n";
    }
}

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene (@genes) { 
    my @treefams   = sort keys %{ $data_ref->{'gene'}->{$gene}->{'treefam'} };
    my $tfam_count = @treefams;
    if ( $tfam_count == 1 ) { 
        print "$gene\t$treefams[0]\n";
    }
    if ( $tfam_count > 1 ) { 
        my $treefam_string = join ', ', @treefams;
        warn "Inconsistent mapping: $gene to $treefam_string\n";
    }
    if ( $tfam_count < 1 ) {
        die "Somehow, gene $gene does not have even one TreeFam family mapped to it!\n";
    }
}

