#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::Util qw(min max);

while ( my $input = <> ) {
    chomp $input;

    my $gene     = q{};
    my $aroian   = 0;
    my $baylor   = 0;
    my $keiser   = 0;
    my $oita     = 0;
    my $obscurus = 0;
    my $ratio    = 0;

    if ( $input =~ /\A (\S+) \t /xms ) {
        $gene = $1;
    }
    else {
        die "Can't parse gene name in: $input\n";
    }

    if ( $input =~ / Aroian \s+ \( (\d+) \s+ g\. \) /xms ) {
        $aroian = $1;
    }
    else {
        die "Can't parse Aroian gene count in: $input\n";
    }

    if ( $input =~ / Baylor \s+ \( (\d+) \s+ g\. \) /xms ) {
        $baylor = $1;
    }
    else {
        die "Can't parse Oita gene count in: $input\n";
    }

    if ( $input =~ / Keiser \s+ \( (\d+) \s+ g\. \) /xms ) {
        $keiser = $1;
    }
    else {
        die "Can't parse obscurus gene count in: $input\n";
    }

    if ( $input =~ / Oita \s+ \( (\d+) \s+ g \. \) /xms ) {
        $oita = $1;
    }
    else {
        die "Can't parse Baylor gene count in: $input\n";
    }

    if ( $input =~ / obscurus \s+ \( (\d+) \s+ g\. \) /xms ) {
        $obscurus = $1;
    }
    else {
        die "Can't parse Keiser gene count in: $input\n";
    }

    my @anhuis = ($baylor, $keiser);
    my $anhui_min = min(@anhuis);

    my @non_anhuis = ($aroian, $oita, $obscurus);
    my $non_anhui_max = max(@non_anhuis); 

    $ratio    = ($anhui_min/$non_anhui_max);
    $ratio    = sprintf("%.3f", $ratio);

    print "$gene\t$ratio\n";
}
