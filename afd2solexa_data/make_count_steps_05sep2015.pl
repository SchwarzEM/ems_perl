#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A \S+ \/ (\S+) \z/xms ) { 
        my $readname = $1;
        $readname =~ s/\.filt\.fastq\z//;
        $readname =~ s/\.fastq\z//;
        print "cat $input | count_simple_fastq_residues.pl 1>$readname.count.txt 2>$readname.count.err ;\n";
    }
}


