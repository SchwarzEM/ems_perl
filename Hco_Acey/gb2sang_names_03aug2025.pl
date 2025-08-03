#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while ( my $input = <>) {
    if ( $input =~ /\A [>] /xms ) {
        if ( $input =~ /\A [>] (\S+) .+ [:][ ] (\S+) /xms ) {
            my $genbank = $1;
            my $sanger  = $2;
            $sanger =~ s/[,]\z//;
            print "$genbank\t$sanger\n";
        }
        else {
            die "Cannot parse header: $input\n";
        }
    }
}
