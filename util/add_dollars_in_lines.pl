#!/usr/bin/env perl

# add_dollars_in_lines.pl -- Erich Schwarz <emsch@its.caltech.edu>, 7/21/2011.
# Purpose: simple hack for getting reliable sum of costs in an ASCII text cost list.

use strict;
use warnings;

my $total = 0;
my $cost  = 0;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ / \A [^\$]* [\$] (\d+\.\d+) [^\$]* \z /xms ) { 
        $cost = $1;
        $total += $cost;
    }
    elsif ( $input =~ / [\$] /xms ) {
        die "Can't parse input: $input\n";
    }
}

print "\n", 'Total: $', "$total\n\n";

