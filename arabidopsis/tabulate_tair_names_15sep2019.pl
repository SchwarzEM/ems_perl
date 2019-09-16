#!/usr/bin/env perl

# tabulate_tair_names_15sep2019.pl  -- Erich Schwarz <ems394@cornell.edu>, 9/15/2019.
# Given output from https://www.arabidopsis.org/tools/bulk/genes/index.jsp, make a clean table of TAIR gene IDs, names, and descriptions.

use strict;
use warnings;
use autodie;

my $data_ref;

my $header = "Gene\tPrimary_name\tAll_names\tDescription";

while (my $input = <>) {
    chomp $input;
    # Data columns:
    # Locus Identifier        Gene Model Name Gene Model Description  Gene Model Type Primary Gene Symbol     All Gene Symbols
    # Sample desired input:
    # AT1G09710       AT1G09710.2     Paralog of DRMY1 with unknown function. protein_coding  DRMY1 PARALOG 1 (DP1)   DRMY1 PARALOG 1 (DP1)
    if ( $input =~ /\A (AT (?: 1|2|3|4|5|C|M)  G\d+) \t \S+ \t ([^\t]*) \t \S+ \t ([^\t]*) \t (.+) \z/xms ) { 
        my $gene_id = $1;
        my $desc    = $2;
        my $name    = $3;
        my $othr    = $4;

        $desc =~ s/\A\s+//;
        $desc =~ s/\s+\z//;

        $name =~ s/\A\s+//;
        $name =~ s/\s+\z//;

        $othr =~ s/\A\s+//;
        $othr =~ s/\s+\z//;
        $othr =~ s/\t/; /g;
        $othr =~ s/;(\S+)/; $1/g; 

        print "$header\n" if $header;
        $header = q{};

        print "$gene_id\t$name\t$othr\t$desc\n";
    }
    elsif ( $input =~ /\A (AT\d+G\d+) \t /xms ) { 
        die "Can't parse input: $input\n";
    }
}

