#!/usr/bin/env perl

# make_fq_readsplit.pl -- Erich Schwarz <emsch@caltech.edu>, 11/6/2012.
# Purpose: given directory full of *.pe.fq files, make script to split them into *.R1.se.fq and *.R2.se.fq files.

use strict;
use warnings;

print '#!/bin/bash', "\n\n";

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \.pe\.fq \z /xms ) { 
        my $filestem = $1;
        print "    split_paired_fastq.or.a.pl -q --r1 \"#0\/1\"  --r2 \"#0\/2\" -i $input --o1 $filestem.R1.se.fq --o2 $filestem.R2.se.fq ;\n";
    }
    else { 
        die "Can't parse input: $input\n"
    }
}

print "\n";

