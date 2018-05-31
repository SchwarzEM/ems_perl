#!/usr/bin/env perl

# part_reformat_canonical_gtf.pl -- Erich Schwarz <ems394@cornell.edu>, 5/29/2018.
# Purpose: extract a subset of WormBase 'canonical genes' in GTF format and make it into simple-minded GFF3.

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+ \t [^\t]* \t) (?:mRNA) (\t .+ \t) ID=\S+?;Parent=\S+?;transcript_id=(\S+?);gene_id=(\S+?) \z/xms ) {
        my $most_data1 = $1;
        my $most_data2 = $2;
        my $tx         = $3;
        my $gene       = $4;
        print "$most_data1", "transcript", "$most_data2", 'ID=', "$tx", ';Parent=', "$gene\n";
    }
    elsif ( $input =~ /\A (\S+ \t [^\t]* \t (?:CDS) \t .+ \t) Parent=\S+?;gene_id=\S+?;transcript_id=(\S+?) \z/xms ) {
       	my $most_data =	$1;
       	my $tx 	 = $2;
        print "$most_data", 'ID=', "$tx.cds", ';Parent=', "$tx\n";
    }
}

