#!/usr/bin/env perl

use strict;
use warnings;
use autodie; 

while (my $input = <>) { 
    chomp $input;
    $input =~ s/N. americanus\t/N. americanus\tNam_/;
    $input =~ s/C. elegans\t/C. elegans\tCele_/;
    $input =~ s/C. briggsae\t/C. briggsae\tCbri_/;
    print "$input\n";
}

