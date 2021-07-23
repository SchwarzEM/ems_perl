#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \/ \S+ (_\d+\.fastq)\.gz \z/xms ) { 
        my $type   = $1;
        my $suffix = $2;
        my $output = $type . $suffix;
        print "gunzip $input ;\n";
        $input =~ s/\.gz\z//;
        print "cat $input >> filt0_RNAseq_reads/$output ;\n";
    }
    else {
        die "Cannot parse: $input\n";
    }
}

