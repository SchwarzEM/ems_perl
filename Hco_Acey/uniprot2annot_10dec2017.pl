#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my $header = "Protein_accession\tGene\tSpecies\tProtein_description";

while (my $input = <>) {
    chomp $input;

    if ( $input =~ /\A [>] tr\| (\S+) \| \S+ \s+ (.+) [ ] OS=(.+) [ ] GN=(\S+) /xms ) {
        my $accession   = $1;
        my $description = $2;
        my $species     = $3;
        my $gene        = $4;

        my $annotation  = "$gene\t$species\t$description";

        if ( exists $data_ref->{'accession'}->{$accession}->{'annotation'} ) {
            my $old_annotation = $data_ref->{'accession'}->{$accession}->{'annotation'};
            die "Redundant annotations for $accession: $old_annotation and $annotation\n";
        }
        else {
            $data_ref->{'accession'}->{$accession}->{'annotation'} = $annotation;
        }
    }
    elsif ( $input =~ /\A [>] tr\| (\S+) \| \S+ \s+ (.+) [ ] OS=(.+) /xms ) {
        my $accession   = $1;
        my $description = $2;
        my $species     = $3;
        my $gene        = "unknown:" . $accession;
        
        # trim off extra junk from $species, if it's there; e.g., trailing "OS=(.+) [ ] GN=(\S+)"
        if ( $species =~ /\A (.+?) \s+ \S+ [=] /xms ) {
            $species = $1;
        }

        my $annotation  = "$gene\t$species\t$description";
        
        if ( exists $data_ref->{'accession'}->{$accession}->{'annotation'} ) {
            my $old_annotation = $data_ref->{'accession'}->{$accession}->{'annotation'};
            die "Redundant annotations for $accession: $old_annotation and $annotation\n";
        }
        else {
            $data_ref->{'accession'}->{$accession}->{'annotation'} = $annotation;
        }
    }

    elsif ( $input =~ /\A [>] /xms ) {
        die "Cannot format UniProt-like FASTA header: $input\n";
    }
}

my @accessions = sort keys %{ $data_ref->{'accession'} };
foreach my $accession (@accessions) {
    print "$header\n" if $header;
    $header = q{};
    my $annotation = $data_ref->{'accession'}->{$accession}->{'annotation'};
    my $output = "$accession\t$annotation";
    print "$output\n";
}

