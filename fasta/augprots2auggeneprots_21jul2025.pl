#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while ( my $input = <> ) {
    chomp $input;
    if ( ( $input !~ /\A [>]/xms ) and ( $input =~ /\S/xms ) ) {
        print "$input\n";
    }
    elsif ( $input =~ /\A ([>]\S+)\.t\d+ .* \z/xms ) {
        my $revised = $1;
        print "$revised\n";
    }
    elsif ( $input =~ /\A [>] /xms ) {
        die "Cannot parse input: $input\n";
    }
}

