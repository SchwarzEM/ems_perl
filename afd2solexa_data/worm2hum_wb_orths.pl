#!/usr/bin/env perl

# worm2hum_wb_orths.pl -- Erich Schwarz <emsch@caltech.edu>, 4/28/2012.
# Purpose: get WBGene/ENSP table of worm-human orthologs from 

use strict;
use warnings;

my $wbgene  = q{};
my $humprot = q{};

my $data_ref;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (WBGene\d+) /xms ) { 
        $wbgene = $1;
    }
    elsif ( $wbgene and ( $input =~ /\A Homo [ ] sapiens \s+ ENSEMBL: (ENSP\d+) /xms ) ) { 
        $humprot = $1;
        $data_ref->{'wbgene'}->{$wbgene}->{'humprot'}->{$humprot} = 1;
    }
}

foreach my $wbg (sort keys %{ $data_ref->{'wbgene'} } ) { 
    if ( exists $data_ref->{'wbgene'}->{$wbg}->{'humprot'} ) { 
        foreach my $hmprt ( sort keys %{ $data_ref->{'wbgene'}->{$wbg}->{'humprot'} } ) { 
            print "$wbg\t$hmprt\n";
        }
    }
}

