#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $canonical = q{};
my $subset    = q{};

my %seen = ();

$canonical = $ARGV[0] if $ARGV[0]; 
$subset    = $ARGV[1] if $ARGV[1];

if ( (! $canonical ) or (! $subset ) ) {
    die "Format: subset_from_canonical_list.pl [canonical list] [subset] > subset_from_canonical_list\n";
}

open my $SUBSET, '<', $subset;
while (my $input = <$SUBSET>) {
    chomp $input;
    $seen{$input} = 1;
}
close $SUBSET;

open my $CANONICAL, '<', $canonical;
while (my $input = <$CANONICAL>) {
    chomp $input;
    if ( exists $seen{$input} ) {
        print "$input\n";
    }
}
close $CANONICAL;

