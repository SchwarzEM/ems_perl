#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $total = 0;

my @times = ();

while (my $input = <>) {
    if ( $input =~ / resources_used \. cput [=] (\d+) [:] (\d+) [:] /xms ) {
        my $hours   = $1;
        my $minutes = $2;
        if ( $hours eq '00' ) { 
            $hours = 0;
        }
        if ( $minutes eq '00' ) {
            $minutes = 0;
        }
        my $time = $hours + ( $minutes / 60 );
        push @times, $time;
        $total = $total + $time;
    }
}

@times = sort { $a <=> $b } @times;

print "Total: $total\n";

foreach my $time (@times) {
    print "$time\n";
}
