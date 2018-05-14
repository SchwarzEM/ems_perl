#!/usr/bin/env perl

# get_geometric_mean.pl -- Erich Schwarz <ems394@cornell.edu>, 11/16/2013.
# Purpose: get a geometric mean of a list of numbers.

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
use Statistics::Descriptive;

my @series = ();

while (my $input = <>) {
    chomp $input;
    if (! looks_like_number($input)) {
        die "Non-numerical value: $input\n";
    }
    push @series, $input;
}

my $stat1 = Statistics::Descriptive::Full->new();

$stat1->add_data(@series);

my $geometric_mean = $stat1->geometric_mean();

print "$geometric_mean\n";

