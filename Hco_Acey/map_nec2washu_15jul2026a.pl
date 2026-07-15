#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while ( my $input = <> ) {
    chomp $input;
    if ( $input =~ /\A (?: [^\t]* \t){9} [^\t]* Parent [=] (\S+) \.t\d+ [;] [^\t]* \t (?: [^\t]* \t){9} [^\t]* transcript [:] (\S+) \z/xms ) {
        my $nec   = $1;
        my $washu = $2;
        print "$nec\t$washu\n";
    }
    else {
        die "Cannot parse: $input\n";
    }
}
