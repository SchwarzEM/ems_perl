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
        print "$ogroup",
              "\t", 
              'min(Aroian,Ilik2,Oita)/max(Anhui,Keiser)',
              "\t",
              'min(Anhui,Keiser)/max(Aroian,Ilik2,Oita)',
              "\t",
              'min(Aroian,Ilik2,Oita,obscurus)/max(Anhui,Keiser)',
              "\t",
              'min(Anhui,Keiser)/max(Aroian,Ilik2,Oita,obscurus)',
              "\t",
              "$annots",
              "\n"
              ;
    }
    else {
        my $necanhui1 = $vals[2];
        my $necanhui  = $vals[3];
        my $necaroian = $vals[4];
        my $neccar1   = $vals[6];
        my $necjapan  = $vals[7];
        my $neckeiser = $vals[8];
        my $nectype3  = $vals[10];
        my @anhuis    = ($necanhui1, $necanhui, $neckeiser);
        my @americs   = ($necaroian, $neccar1, $necjapan);
        my @necators  = ($necaroian, $neccar1, $necjapan, $nectype3);

        my $anhui_max = max(@anhuis);
        my $anhui_min = min(@anhuis);

        my $americ_max = max(@americs);
        my $americ_min = min(@americs);

        my $necator_max = max(@necators);
        my $necator_min = min(@necators);

        my $ratio1 = ( $americ_min / ( $anhui_max  + 0.001 ) );
        my $ratio2 = ( $anhui_min  / ( $americ_max + 0.001 ) );

        my $ratio3 = ( $necator_min / ( $anhui_max   + 0.001 ) );
        my $ratio4 = ( $anhui_min   / ( $necator_max + 0.001 ) );

        print "$ogroup\t$ratio1\t$ratio2\t$ratio3\t$ratio4\t$annots\n";
    }
}

