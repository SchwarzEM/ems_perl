#!/usr/bin/env perl

use strict;
use warnings;

print '#!/bin/bash', "\n\n";

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A Gene ( (?: \t \S+)+ ) \z/xms ) { 
        my @fields = split /\t/, $input;

        my $slice_number = @fields;

        # Note: the number which we want to give 'cut' is 1 more than the number we want to use to pull text out of @fields.
        foreach my $i (2..$slice_number) { 
            # Example of a type: "log10(24.PI/L3i)"
            my $output_type = $fields[ ($i-1) ];
            $output_type    =~ s/\//.vs./g;
            $output_type    =~ s/\(/./g;
            $output_type    =~ s/\)/./g;
            $output_type    =~ s/[.]+/./g;
            $output_type    =~ s/\A[.]//;
            $output_type    =~ s/[.]\z//;
            print "    cut -f 1,$i Acey_2012.10.24.pgene_log10_X.vs.Y_24mar2013.txt > Acey_2012.10.24.pgene_", $output_type, "_24mar2013.txt ;\n";
        }
    }
}

print "\n";

