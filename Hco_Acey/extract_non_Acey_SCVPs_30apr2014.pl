#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = "Name\tSpecies\tDescription";

while (my $input = <>) { 
    chomp $input;

    print "$header\n" if $header;
    $header = q{};

    if ( $input =~ /\A > (Cbri_\S+) .+ (WBGene\d+) /xms ) { 
        my $name = $1;
        my $desc = $2;
        print "$name\tCaenorhabditis briggsae\t$desc\n";
    }
    elsif ( $input =~ /\A > (Cele_\S+) .+ (WBGene\d+) /xms ) {
        my $name = $1;
        my $desc = $2;
        print "$name\tCaenorhabditis elegans\t$desc\n";
    }
    elsif ( $input =~ /\A > (Hbac_\S+) .+ \s (\S+) \s* \z/xms ) {
        my $name = $1;
        my $desc = $2;
        print "$name\tHeterorhabditis bacteriophora\t$desc\n";
    }
    elsif ( $input =~ /\A > (Hcon_\S+) .+ \s (\S+) \s* \z/xms ) {
        my $name = $1;
        my $desc = $2;
        print "$name\tHaemonchus contortus\t$desc\n";
    }
    elsif ( $input =~ /\A > (Name_\S+) .+ \s (\S+) \s* \z/xms ) {
        my $name = $1;
        my $desc = $2;
        print "$name\tNecator americanus\t$desc\n";
    }
    elsif ( $input =~ /\A > (Ppac_\S+) .+ (WBGene\d+) /xms ) {
        my $name = $1;
        my $desc = $2;
        print "$name\tPristionchus pacificus\t$desc\n";
    }
}

