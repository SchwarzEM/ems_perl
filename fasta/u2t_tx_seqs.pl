#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while ( my $input = <> ) {
    chomp $input;
    if ( $input !~ /\A [>]/xms ) {
        $input =~ s/U/T/g;
    }
    print "$input\n";
}
