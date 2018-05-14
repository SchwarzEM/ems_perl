#!/usr/bin/env perl

# write_blankace.pl -- Erich Schwarz <emsch@its.caltech.edu>, 6/24/2008.
# Purpose: given a .ace file, write a simple .ace file wiping out previous data.

use strict;
use warnings;

my %write_outs = ();

while (my $input = <>) { 
    if ( $input =~ / \A 
                     ( \S+ \s : \s 
                       \" [^\"]+ \" ) /xms ) { 
        my $entry = $1;
        $write_outs{$entry} = 1;
    }
}

print "\n";

foreach my $entry (sort keys %write_outs) { 
    print "-D $entry\n";
    print "\n";
}

