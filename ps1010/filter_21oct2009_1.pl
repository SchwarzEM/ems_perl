#!/usr/bin/env perl

# filter_21oct2009_1.pl -- filter only "Transcript \"[^\"\s]+\" lines from huge GFF files.

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ / Transcript \s+ \" [^\"\s]+ \" /xms ) { 
        print "$input\n";
    }
} 

