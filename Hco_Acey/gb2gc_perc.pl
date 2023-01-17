#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $locus = q{};

my $fullA = 0;
my $fullC = 0;
my $fullG = 0;
my $fullT = 0;

while (my $input = <>) {
    chomp $input;

    # Sample input:
    # LOCUS       Necator_chrI_37100-46090   8991 bp  DNA
    if ( $input =~ /\A LOCUS \s+ (\S+) \s+ \d+ \s+ bp \s+ DNA \s* \z/xms ) {
        $locus = $1;
    }
    # Sample input:
    # Sample input:   BASE COUNT     2896 a   1909 c  1771 g   2415 t
    # Another sample: BASE COUNT     2296 a   1511 c  1411 g   2113 t   1035 n
    elsif ( $input =~ /\A BASE \s COUNT \s+ (\d+) \s+ a \s+ (\d+) \s+ c \s+ (\d+) \s+ g \s+ (\d+) \s+ t \b /xms ) {  
        my $A = $1;
        my $C = $2;
        my $G = $3;
        my $T = $4;

        $fullA = $fullA + $A;
        $fullC = $fullC + $C;
        $fullG = $fullG + $G;
        $fullT = $fullT + $T;

        my $all = ( $A + $C + $G + $T );
        my $GC  = ( $C + $G );
        my $GC_perc = 100 * ( $GC / $all );

        print "$locus\t$GC_perc\n";
        $locus = q{};
    }
    else {
        die "Cannot parse: \"$input\"\n";
    }
}

my $full_all = ( $fullA + $fullC + $fullG + $fullT );
my $full_GC  = ( $fullC + $fullG );
my $full_GC_perc = 100 * ( $full_GC / $full_all );

print "LOCUS_total\t$full_GC_perc\n";

