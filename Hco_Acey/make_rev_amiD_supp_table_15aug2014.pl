#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $pgrp_metazoa  = $ARGV[0];
my $orig_table    = $ARGV[1];
my $header        = "Phylogenetic sequence name\tCategory\tSpecies\tProtein accession no.\tGene name\tNon-bacterial phylogeny";
my $nonbact_phylo = q{};
my %pgrp_seen     = ();

open my $PGRP, '<', $pgrp_metazoa;
while (my $input = <$PGRP>) {
    chomp $input;
    while ( $input =~ / ['] ([^'\s\/]+) \/ \d+ [-] \d+ ['] /xmsg ) {
        my $seen_seq = $1;
        $pgrp_seen{$seen_seq} = 1;
    }
}
close $PGRP;

open my $ORIG, '<', $orig_table;
while (my $input = <$ORIG>) {
    chomp $input;
    if ( $input =~ /\A ([^\t]+) \t ([^\t]+) \t ([^\t]* \t [^\t]* \t [^\t]*) \z/xms ) { 
        my $phylo_name = $1;
        my $category   = $2;
        my $other_text = $3;
        $nonbact_phylo = q{};

        # Fix this inconsistency:
        if ( $category eq 'Metagenome') {
            $category = 'Metagenomes';
        }

        if ( $category ne 'Bacteria' ) { 
            if ( exists $pgrp_seen{$phylo_name} ) {
                $nonbact_phylo = 'PGRP';
            }
            elsif ( ( $category ne 'Metagenomes' ) and ( $category ne 'Viruses' ) ) {
                $nonbact_phylo = 'HGT';
            }
        }
        print "$header\n" if $header;
        $header = q{};
        print "$phylo_name\t$category\t$other_text\t$nonbact_phylo\n";
    }
}

