#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t \. \t (.+) \z/xms ) {
        my $front = $1;
        my $back  = $2;
        print "$front\tnot_lifted_over\t$back\n";
    }
    elsif ( $input =~ /\A (\S+) \t \S+ \t (.+) \z/xms ) { 
        my $front = $1;
        my $back  = $2;
        print "$front\tliftOver-UCSC+N2_WS264 \t$back\n";
    }
    elsif ( $input !~ /\A [#] /xms ) {
       die "Cannot parse: $input\n";
    }
}

