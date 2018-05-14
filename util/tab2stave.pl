#!/usr/bin/env perl

# tab2stave.pl -- Erich Schwarz <emsch@its.caltech.edu>, 4/25/2008.

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    $input =~ s/\t/|/g;
    print "$input\n";
}

