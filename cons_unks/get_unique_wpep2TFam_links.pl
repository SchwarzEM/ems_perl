#!/usr/bin/env perl

# get_unique_wpep2TFam_links.pl -- Erich Schwarz <ems394@cornell.edu>, 8/14/2013.
# Purpose: given tabular psi-BLAST of largest wormpep isoforms (i.e., w/ 1 protein-gene) versus TreeFam etc., find all *unique* links (ignoring multiple ones).

use strict;
use warnings;

my $data_ref;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \s+ (TF\d+) /xms ) { 
        my $wpep = $1;
        my $tfam = $2;
        $data_ref->{'wpep'}->{$wpep}->{'tfam'}->{$tfam} = 1;
    }
}

my @final_wpeps = sort keys %{ $data_ref->{'wpep'} };
foreach my $wpep (@final_wpeps) {
    my @final_tfams = sort keys %{ $data_ref->{'wpep'}->{$wpep}->{'tfam'} };
    my $tfam_count  = @final_tfams;
    if ( $tfam_count == 1 ) {
        print "$wpep\t$final_tfams[0]\n";
    }
}

