#!/usr/bin/env perl

# wormpep2wbgene_table.pl -- Erich Schwarz <ems394@cornell.edu>, 6/24/2018.
# Purpose: given headers from wormpep (in format seen for WS264, 4/2018) emit "gene\tCDS\tlocus" TSV table; do not emit 2+ identical IDs.

use strict;
use warnings;

my %seen         = ();

# Typical wormpep190 lines:
# >4R79.1b        CE39659 WBGene00003525  locus:nas-6 ..
# >AC7.3  CE07653 WBGene00014997

while (my $input = <>) { 
    chomp $input;
    if ($input =~ /\A [>] (\S+) \b .+ (WBGene\d+) /xms) { 
        my $cds         = $1;
        my $gene        = $2;
        my $locus       = q{};

        $cds =~ s/[a-z]\z//;

        if ( $input =~ / locus = (\S+) /xms ) {
            $locus = $1;
        }

        if (! exists $seen{$gene} ) {
            print "$gene\t$cds\t$locus\n";
            $seen{$gene} = 1;
        }
    }
    elsif ( $input =~ /\A [>] /xms ) { 
        die "Cannot parse FASTA header line: $input\n";
    }
}

