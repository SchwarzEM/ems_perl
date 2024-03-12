#!/usr/bin/env perl

use strict;
use autodie;
use warnings;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A [>] (\S+) (.*) \z/xms ) {
        my $seq  = $1;
        my $comm = $2;
        $seq =~ s/(\.\d+[a-z]*)\.\d+\z/$1/;
        $input = '>' . $seq . $comm;
    }
    print "$input\n";
}
