#!/usr/bin/env perl

# oocalc_csv2tsv.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/27/2010.
# Purpose: given a CSV from oocalc festooned with unwanted \" text delimiters, strip them out.

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    $input =~ s/\A\"//;
    $input =~ s/\"\z//;
    $input =~ s/\t\"/\t/g;
    $input =~ s/\"\t/\t/g;
    print "$input\n";
}


