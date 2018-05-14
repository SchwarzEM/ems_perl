#!/usr/bin/env perl

# strip_wbg_name_alts.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/24/2010.
# Purpose: remove the staved extra terms from WBGene\d+ ID that starts a line of text.

use strict;
use warnings;

my $name = q{};
my $text = q{};

while (my $input = <> ) {
    chomp $input;
    if ( $input =~ / \A (WBGene\d+\|\S+) (.*) \z /xms ) { 
        $name = $1;
        $text = $2;
        $name =~ s/(WBGene\d+)\|\S+/$1/;
        $input = $name . $text;
    }
    print "$input\n";
}


