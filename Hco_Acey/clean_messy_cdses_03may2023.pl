#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (.+) \z/xms ) { 
        my $id   = $1;
        my $text = $2;
        if ( $id ne 'Gene' ) {
            $id =~ s/[A-Za-z]\z//;
        }
        if ( $id =~ /\A (\S+\.\d+)\.\d+ \z/xms ) {
            $id = $1;
        }
        # Some names are: 'F53F10.2b.1'
        if ( $id =~ /\A (\S+\.\d+)[a-z]\.\d+ \z/xms ) {
            $id = $1;
        }
        print "$id\t$text\n";
    }
    else {
        die "Cannot parse: $input\n";
    }
}
