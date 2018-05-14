#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = '#!/bin/bash';

while (my $input = <>) {
    chomp $input;
    # Sample input:
    # ATML1__LGO.vs.lgo-2.upreg.2016.02.06.func_bool_precursor.txt
    if ( $input =~ /\A (\S+)\.func_bool_precursor\.txt \z/xms ) {
        my $stem = $1;
        my $output = "$stem.func_bool.txt";
        print "$header\n\n" if $header;
        $header = q{};
        print "    qual2func_table_arath.pl -t $input -g gene_association.05feb2016.tair";
        print ' 1>';
        print "$output 2>>test.err ;\n";
    }
    else {
        die "Cannot parse input: $input\n";
    }
}

print "\n" if (! $header);

