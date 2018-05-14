#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) ((?: \t [^\t]* ){10}) ((?: \t [^\t]* ){15}) \z/xms ) {
        my $gene = $1;
        my $sig  = $2;
        my $tpms = $3;
        if ( $sig =~ /\S/xms ) {
            print "$gene$tpms\n";
        }
    }
    else {
       die "Cannot parse input line: $input\n";
    }
}

