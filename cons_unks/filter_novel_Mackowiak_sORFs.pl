#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A [#] /xms ) {
        print "$input\n";
    }
    elsif ( $input =~ /\A (?: [^\t]* \t){15} (\S+) \t (?: [^\t]* \t){3} True \t True \t /xms ) {
        my $status = $1;
        if ( $status ne 'annotated' ) { 
            print "$input\n";
        }
    }
}

