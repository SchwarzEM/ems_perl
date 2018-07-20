#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \s+ \S+ \s+ 
                       (\S+) \s+ 
                       (\S+) \s+ \S+ \s+ \d+ \s+ \d+ \s+ 
                       (\d+) \s+ 
                       (\d+) \s+ 
                       (\S+) \s+ \S+ \s+ \S+ \s+ \S+ \s+ \S+ \s+ \S+ \s+ 
                       (\S+) /xms 
       ) { 
        my $seq        = $1;
        my $rfam_name  = $2;
        my $rfam_id    = $3;
        my $orig_start = $4;
        my $orig_end   = $5;
        my $strand     = $6; 
        my $e_val      = $7;

        my @coords = ($orig_start, $orig_end);
        @coords = sort { $a <=> $b } @coords;
        my $start_nt = $coords[0];
        my $end_nt   = $coords[1];
        $start_nt--;

        print "$seq";
        print "\t";
        print "$start_nt";
        print "\t";
        print "$end_nt";
        print "\t";
        print "$rfam_name";
        print "_";
        print "$rfam_id";
        print "\t";
        print "$e_val";
        print "\t";
        print "$strand";
        print "\n";
    }
    else {
        die "Cannot parse: $input\n";
    }
}


