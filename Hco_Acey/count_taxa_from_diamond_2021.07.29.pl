#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $infile   = q{};
$infile      = $ARGV[0] if $ARGV[0];
my %taxcount = ();

open my $INFILE, '<', $infile;
while (my $input = <$INFILE>) {
    chomp $input;
    if ( $input =~ /OS=(\S.+\S)\s+OX=/xms ) {
        my $taxon = $1;
            $taxcount{$taxon}++;
    }
    else {
        warn "Cannot parse taxon in: $input\n";
    }
}
close $INFILE;

my @taxa = sort keys %taxcount;

foreach my $taxon (@taxa) {
    my $count = $taxcount{$taxon};
    print "$count\t$taxon\n";
}

