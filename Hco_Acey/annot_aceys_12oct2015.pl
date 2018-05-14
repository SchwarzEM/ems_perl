#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %seq2annot = (
    'Acey_s0184.g991.t1'  => "ASP-2 [positive control]",
    'Acey_s0045.g1091.t3' => "MTP-1 [positive control]",
    'Acey_s0045.g1116.t1' => "paralog of MTP-1 [positive control]",
    'Acey_s0010.g910.t3'  => "Mannose receptor [highest priority]",
    'Acey_s0004.g1962.t1' => "Asialoglycoprotein receptor",
    'Acey_s0230.g2988.t1' => "Mannose receptor [lower priority]",
    'Acey_s0001.g225.t1'  => "Zinc protease, upregulated during infection",
);

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A > (\S+) /xms ) { 
        my $seq = $1;
        if ( exists $seq2annot{$seq} ) {
            $input =~ s/\A>$seq/>$seq  $seq2annot{$seq}  /;
        }
    }
    print "$input\n";
}

