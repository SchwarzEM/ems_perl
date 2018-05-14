#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    my $briggsae_count     = 0;
    my $briggsae_alt_count = 0;

    while ( $input =~ /\(briggsae\)/xmsg ) { 
        $briggsae_count++;
    }
    while ( $input =~ /\(briggsae_alt\)/xmsg ) {
        $briggsae_alt_count++;
    }

    if ( ( $briggsae_count == 1 ) and ( $briggsae_alt_count == 1 ) ) {
        print "$input\n";
    }
}
