#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = '#!/bin/bash';

while (my $input = <>) {
    chomp $input;
    # Sample input:
    # ATML1__LGO_atml1-3.vs.atml1-3.upreg.gene_list.txt
    if ( $input =~ /\A (\S+)_list\.txt \z/xms ) {
        my $stem = $1;
        $stem =~ s/\.gene\z//;
        my $output = "$stem.2016.03.02.func_bool_precursor.txt";
        print "$header\n\n" if $header;
        $header = q{};
        print "    build_func_boolean.pl -t -p $input -g TAIR10_cds_20101214_updated.gene_list.txt";
        print ' 1>';
        print "$output 2>>test.err ;\n";
    }
    else {
        die "Cannot parse input: $input\n";
    }
}

print "\n" if (! $header);
