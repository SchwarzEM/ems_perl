#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::Util qw(min max);

while ( my $input = <> ) {
    chomp $input;
    my @vals = split '\t', $input;
    my $ogroup  = q{};
    my $annots  = q{};
    my $ranking = 'NA';

    if ( $input =~ /\A (\S+) \t (.+) \z/xms ) {
        $ogroup = $1;
        $annots = $2;
    }
    else {
        die "Cannot parse: $input\n";
    }

    if ( $ogroup eq 'Orthogroup' ) {
        print "$ogroup\tAroian.Ilik2.Oita/Anhui,Keiser\t$annots\n";
    }
    else {
        my $necanhui1 = $vals[2];
        my $necanhui  = $vals[3];
        my $necaroian = $vals[4];
        my $neccar1   = $vals[6];
        my $necjapan  = $vals[7];
        my $neckeiser = $vals[8];
        my @anhuis     = ($necanhui1, $necanhui, $neckeiser);
        my @non_anhuis = ($necaroian, $neccar1, $necjapan);

        my $anhui_val     = max(@anhuis);
        my $non_anhui_val = min(@non_anhuis);
        $anhui_val = ($anhui_val + 0.001);
        my $ratio = ( $non_anhui_val / $anhui_val );
        print "$ogroup\t$ratio\t$annots\n";
    }
}

