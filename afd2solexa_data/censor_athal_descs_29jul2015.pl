#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) { 
    chomp $input;
    $input =~ s/BEST \s Arabidopsis \s thaliana \s protein \s match \s is: [^;]+ ;//xms;
    $input =~ s/Has \s \d+ \s Blast \s hits \s to \s \d+ \s proteins \s in \s \d+ \s species: .+ \(source: \s NCBI \s BLink\)//xms;
    $input =~ s/[ ]+/ /g;
    $input =~ s/; ;/;/g;
    $input =~ s/; \././g;
    print "$input\n";
}

