#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $peppers = q{};
my $genes   = q{};

$peppers = $ARGV[0] if $ARGV[0];
$genes   = $ARGV[1] if $ARGV[1];

my %seen = ();

if ( (! $peppers ) or (! $genes ) ) {
    die "Format: tab_peppers_04jun2021.pl [pepper genes] [full gene list] > [pepper table]\n"; 
}

open my $PEPPER, '<', $peppers;
while (my $input = <$PEPPER>) {
    chomp $input;
    $seen{$input} = 1;
}
close $PEPPER;

open my $GENES, '<', $genes;
while (my $input = <$GENES>) {
    chomp $input;
    if ( exists $seen{$input} ) {
        print "$input\tPepper_gene\n";
    }
    else {
        print "$input\tNon_pepper_gene\n";
    }
}
close $GENES;
