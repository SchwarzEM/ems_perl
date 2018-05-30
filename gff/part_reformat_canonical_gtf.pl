#!/usr/bin/env perl

# part_reformat_canonical_gtf.pl -- Erich Schwarz <ems394@cornell.edu>, 5/29/2018.
# Purpose: extract a subset of WormBase 'canonical genes' in GTF format and make it into simple-minded GFF3.

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+ \t [^\t]* \t (?:transcript) \t .+ \t) gene_id\s\"(\S+?)\";\stranscript_id\s\"(\S+?)\";/xms ) {
        my $most_data = $1;
        my $gene      = $2;
        my $tx        = $3;
        print "$most_data", 'ID=', "$tx", ';Parent=', "$gene\n";
    }
    elsif ( $input =~ /\A (\S+ \t [^\t]* \t (?:CDS) \t .+ \t) gene_id\s\"\S+?\";\stranscript_id\s\"(\S+?)\";/xms ) {
       	my $most_data =	$1;
       	my $tx 	 = $2;
        print "$most_data", 'ID=', "$tx.cds", ';Parent=', "$tx\n";
    }
}
