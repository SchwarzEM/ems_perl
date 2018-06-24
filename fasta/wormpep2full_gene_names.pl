#!/usr/bin/env perl

# wormpep2full_gene_names.pl -- Erich Schwarz <ems394@cornell.edu>, 4/3/2018
# Purpose: given headers from wormpep (in format seen for WS264, 4/2018) emit full ID stream; do not emit 2+ identical IDs.

use strict;
use warnings;

my %cds2fullname = ();
my %possible_cds = ();
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

        my @id_tags = ($gene, $cds);

        if ( $input =~ / locus = (\S+) /xms ) {
            $locus = $1;
            push @id_tags, $locus;
        }

        my $full_id = join '|', @id_tags;

        if (! exists $seen{$full_id} ) {
            print "$full_id\n";
            $seen{$full_id} = 1;
        }
    }
    elsif ( $input =~ /\A [>] /xms ) { 
        die "Cannot parse FASTA header line: $input\n";
    }
}

