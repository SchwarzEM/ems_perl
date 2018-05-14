#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = "Name\tSpecies\tDescription\tAccession";

while (my $input = <>) { 
    chomp $input;
    # Sample input:
    # >Tadh_XP_002118895.1  gi|196018934|ref|XP_002118895.1| hypothetical protein TRIADDRAFT_62885 [Trichoplax adhaerens]
    if ( $input =~ /\A > (\S+) \s+ (\S+) \s+ (\S.+\S) \s+ \[ (\S.+\S) \] \s* \z/xms ) { 
        my $name    = $1;
        my $acc     = $2;
        my $desc    = $3;
        my $species = $4;
        print "$header\n" if $header;
        $header = q{};
        print "$name\t$species\t$desc\t$acc\n";
    }
    else { 
        die "Can't parse input: $input\n";
    }
}

