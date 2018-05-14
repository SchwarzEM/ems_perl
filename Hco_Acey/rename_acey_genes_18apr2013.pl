#!/usr/bin/env perl

# rename_acey_genes_18apr2013.pl -- Erich Schwarz <ems394@cornell.edu>, 4/18/2013.
# Purpose: systematically rename sequences / genes / transcripts / etc. from Acey_2012.08.05_X to Acey_sX (to make names easier to say or read).

use strict;
use warnings;

while (my $input = <>) { 
    $input =~ s/Acey_2012\.08\.05_/Acey_s/g;
    print $input;
}

