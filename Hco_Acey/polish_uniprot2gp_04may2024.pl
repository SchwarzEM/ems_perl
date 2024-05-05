#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (.* QR680_\d+ .*) /xms ) {
        my $uniprot = $1;
        my $gb_text = $2;
        while ( $gb_text =~ /(QR680_\d+)/xmsg ) {
            my $gb = $1;
            print "$uniprot\t$gb\n";
        }
    }
}
