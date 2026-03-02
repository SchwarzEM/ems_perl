#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $bed = q{};
$bed    = $ARGV[0] if $ARGV[0];

if (! $bed ) {
    die "Format: label_4th_bed_col_02mar2026.pl [input BED] > [output BED]\n";
}

open my $BED, '<', $bed;
while ( my $input = <$BED> ) { 
    chomp $input;
    # Sample input:
    # Necator_chrIV     76      1444    .       .       -       AUGUSTUS        gene    .       ID=Necator_chrIV.g13292
    if ( $input =~ /\A (\S+ \t \d+ \t \d+ \t) \S+ (\t \S+ \t \S+ \t \S+ \t gene \t \S+ \t ID[=] (\S+)) \z/xms ) {
        my $front = $1;
        my $back  = $2;
        my $gene  = $3;
        my $output = "$front$gene$back";
        print "$output\n";
    }
    else {
       die "From input BED file $bed, cannot parse: $input\n";
    }
}
close $BED;

