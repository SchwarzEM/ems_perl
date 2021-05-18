#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %entry2number = ();

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) (\t .+) \z/xms ) {
        my $entry = $1;
        my $data  = $2;
        $entry2number{$entry}++;
        my $number = $entry2number{$entry};
        print "$entry.$number\t$data\n";
    }
    else {
        die "Cannot parse: $input\n";
    }
}
