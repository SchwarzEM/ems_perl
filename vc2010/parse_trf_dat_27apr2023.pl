#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = "Sequence\tStart_nt\tEnd_nt\tTR_len\tUnit_len\tUnit_seq";

my $sequence = q{};

while (my $input = <> ) {
    chomp $input;
    if ( $input =~ /\ASequence [:] \s+ (\S+)\z/xms ) {
        $sequence = $1;
    }
    # Sample input:
    # 1 476 6 79.7 6 97 0 758 18 18 31 31 1.96 GGCTTA GGCTTAGGC[...]
    elsif ( $input =~ /\A (\d+) \s+ (\d+) \s+ (\d+) \s+ .+ \s+ ([ACGT]+) \s+ [ACGT]+ \s* \z/xms ) {
        my $start_nt = $1;
        my $end_nt   = $2;
        my $period   = $3;
        my $repeat   = $4;
        my $length = ($end_nt - $start_nt) + 1;
        if ( $sequence =~ /\S/xms ) {
            if ( $header ) {
                print "$header\n";
                $header = q{};
            }
            print "$sequence\t$start_nt\t$end_nt\t$length\t$period\t$repeat\n";
        }
        else {
            die "No sequence ID, so cannot print: $input\n";
        }
    }
    elsif ( $input =~ /\A (\d+) \s+ (\d+) \s+ (\d+) /xms ) {
        die "Unparsed input with sequence $sequence: $input\n";
    }
}
