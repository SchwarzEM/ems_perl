#!/usr/bin/env perl

# wpep2table_v01.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/25/2009.
# Purpose: get orderly table of peptide, CDS, and WBGene names from WS200 proteome FASTA headers.

use strict;
use warnings;

my $cds     = q{};
my $wbgene  = q{};
my $prot    = q{};
my $uniprot = q{};
my $cgc     = q{};
my $ncbi    = q{};

my %cds2data = ();

while (my $input = <>) { 
    if ($input =~ / \A > /xms ) { 
        store_data();
        clear_data();
    }
    if ($input =~ / \A > (\S+) \s+ ([A-Z]+\d+) \s+ (WBGene\d+) /xms ) {
        $cds    = $1;
        $prot   = $2;
        $wbgene = $3;
        $cgc    = $cds;
        $cgc    =~ s/[a-z]\z//;
    }
    if ($input =~ / \A > .+ locus:(\S+) /xms ) {
        $cgc = $1;
    }

    if ($input =~ / \A > .+ UniProt:(\S+) /xms ) {
        $uniprot = $1;
    }
    if ($input =~ / \A > .+ protein_id:(\S+) /xms ) {
        $ncbi = $1;
    }
}
store_data();

print "WBGene\tGene\tCDS\tProt.\tUniProt\tNCBI\n";

foreach my $cds1 (sort { $cds2data{$a} cmp $cds2data{$b} } keys %cds2data) { 
    print "$cds2data{$cds1}\n";
}

sub store_data { 
    if ($cds) { 
        $cds2data{$cds} = "$wbgene\t$cgc\t$cds\t$prot\t$uniprot\t$ncbi";
    }
}

sub clear_data { 
    $cds     = q{};
    $wbgene  = q{};
    $prot    = q{};
    $uniprot = q{};
    $cgc     = q{};
    $ncbi    = q{};
}

