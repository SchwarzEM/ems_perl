#!/usr/bin/env perl

# Mihoko_TF_table.pl -- Erich Schwarz <emsch@its.caltech.edu>, 6/5/2011.
# Purpose: make simple table listing past Mihoko TFs.

use strict;
use warnings;

my $gene = q{};

print "Gene\tScreened_in_KS2009\n";

while (my $input = <>) { 
    chomp $input;
    $input =~ s/\s*\z//;
    print "$input", "\t", 'KS2009', "\n", ;
}


