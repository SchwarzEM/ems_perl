#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;
my $sequence = q{};

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > (\S+) /xms ) { 
        $sequence = $1;
    }
    else { 
        $input =~ s/\s//g;
        if ( $sequence and ( $input =~ /\S/ ) ) {
            $data_ref->{'seq'}->{$sequence} .= $input;
        }
    }
}

my @sequences = sort keys %{ $data_ref->{'seq'} };

my $header = "CDS\tMotif_density\tMotifs\tSeq_length\n";

LOOP: foreach my $seq (@sequences) { 
    my $residues = $data_ref->{'seq'}->{$seq};
    my $length   = length $residues;
    my $motifs   = 0;

    if ( $length == 0 ) {
        next LOOP;
    }
    while ($residues =~ /RFG[KR][KR]/g) {
        $motifs++;
    }
    my $density = ($motifs/$length);
    $density    = sprintf("%.5f", $density);

    if ($motifs >= 1) {
        print $header if $header;
        $header = q{};
        print "$seq\t$density\t$motifs\t$length\n";
    }
}

