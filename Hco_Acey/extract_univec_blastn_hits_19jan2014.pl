#!/usr/bin/env perl

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;

    # Sample input, which Perl is having a horrible time parsing properly:
    # qseqid                    sseqid                          pident  length mismatch gapopen qstart  qend    sstart send evalue bitscore
    # Acey_rnaseq_cDNA_000002	gnl|uv|AF327711.1:4497-4713	100.00	17	0	0	3186	3202	130	146	   89	34.2

    my @values   = split /\s+/, $input;
    my $qseqid   = $values[0];
    my $qstart   = $values[6];
    my $qend     = $values[7];
    my $bitscore = $values[11];
    if ( $bitscore >= 19 ) { 
        print "$qseqid\t$qstart\t$qend\n"; 
    }
}

