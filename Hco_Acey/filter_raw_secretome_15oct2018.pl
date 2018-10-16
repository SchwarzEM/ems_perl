#!/usr/bin/env perl

# filter_raw_secretome_15oct2018.pl -- Erich Schwarz <ems394@cornell.edu>, 10/15/2018.
# Purpose: filter and fix up raw results for Acey genes from http://www.cbs.dtu.dk/services/SecretomeP

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input !~ /signal peptide predicted by SignalP/ ) {
        if ( $input =~ /\A (\S+) \s+ (\S+) \s/xms ) {
            my $gene     = $1;
            my $nn_score = $2;
            if ( $nn_score >= 0.6 ) {
                $gene = "Acey_s" . $gene;
                $gene =~ s/\.t\z//;
                $gene =~ s/\.\z//;
                print "$gene\n";
            }
        }
    }
}
