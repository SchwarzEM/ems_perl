#!/usr/bin/env perl

# rev_ens2omim.pl -- Erich Schwarz <emsch@its.caltech.edu>, 4/28/2012.
# Purpose: given naive EnsMart download from ensembl, make a better version which only has disease and ENSP entries.

use strict;
use warnings;

my $prot     = q{};
my $omim     = q{};
my $omim_txt = q{};

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A ENSG\d+ \t ENST\d+ \t ([^\t]+) \t (\d+) \t (ENSP\d+) \z/xms ) { 
        $omim_txt = $1;
        $omim     = $2;
        $prot     = $3;
        print "$prot\t$omim\t$omim_txt\n";
    }
}


