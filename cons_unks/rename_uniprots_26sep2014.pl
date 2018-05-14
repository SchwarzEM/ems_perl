#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %id2sp = (
    3702   => 'arabidopsis',
    6239   => 'elegans',
    7227   => 'drosophila',
    9606   => 'human',
    10090  => 'mouse',
    284812 => 'pombe',
    559292 => 'cerevisiae',
    7955   => 'zebrafish',
    44689  => 'dictyostelium',
);

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A uniprot_26sep2014_ (\d+) \.dat \z/xms ) {
        my $id = $1;
        if ( exists $id2sp{$id} ) {
            print "    mv -i $input uniprot_26sep2014_" . "$id2sp{$id}" . "_" . "$id.dat ;\n";
        }
    }
}

