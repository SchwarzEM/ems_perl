#!/usr/bin/env perl

# decomment_newick.pl -- Erich Schwarz <ems394@cornell.edu>, 3/28/2013:
# Purpose: given a Newick file with obnoxious comments at the node labels, strip them out without otherwise affecting the file.

use strict;
use warnings;

while (my $input = <>) { 
    $input =~ s/\'(\S+)[^\']+\'/$1/g;
    print $input;
}

