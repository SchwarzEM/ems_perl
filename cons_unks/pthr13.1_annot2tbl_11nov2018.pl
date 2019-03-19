#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

# Sample input line:
# CAEEL|WormBase=WBGene00010091|UniProtKB=Q20815		PTHR22947:SF21	MAJOR SPERM PROTEIN	SPERM SPECIFIC FAMILY, CLASS P-RELATED	

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A [^\t]+ UniProtKB= (\S+) \t [^\t]* \t (\S+) \t/xms ) {
        my $uniprot = $1;
        my $panther = $2;
        $data_ref->{'uniprot'}->{$uniprot}->{'panther'}->{$panther} = 1;
    }
    else {
        die "Cannot parse input: $input\n";
    }
}

my @uniprots = sort keys %{ $data_ref->{'uniprot'} };
foreach my $uniprot (@uniprots) {
    my @panthers = sort keys %{ $data_ref->{'uniprot'}->{$uniprot}->{'panther'} };
    my $panther_text = join ';', @panthers;
    print "$uniprot\t$panther_text\n";
}

