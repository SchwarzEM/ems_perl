#!/usr/bin/env perl

use strict;
use warnings;

use List::MoreUtils qw(uniq);

my $data_ref;

my @uniprots = ();
my @pmids    = ();

while (my $input = <>) { 
    chomp $input;
    # Delete trailing space strings.
    $input =~ s/\s+\z//;

    # UniProt uses '//' to denote end of a given record.  Handy!
    if ( $input =~ /\A \/ \/ /xms ) {
        map_print_clear_data();
    }

    # UniProts:
    # AC   P31946; A8K9K2; E1P616;

    elsif ( $input =~ /\A AC \s+ (\S.*\S) /xms ) { 
        my $uniprot_text = $1;
        $uniprot_text =~ s/;//g;
        my @new_uniprots = split /\s+/, $uniprot_text;
        push @uniprots, @new_uniprots;
    }

    # Note that the next elsif is designed to capture PubMed IDs that are not at the *start* of their line.
    # E.g.:  RX   MEDLINE=88125000; PubMed=3257574; DOI=10.1073/pnas.85.4.1184;

    elsif ( $input =~ /\A RX \s .* \s PubMed=(\d+) ; /xms ) { 
        my $pmid = $1;
        push @pmids, $pmid;
    }
}

# Failsafe, but probably not needed:
map_print_clear_data();

sub map_print_clear_data { 
    if (@uniprots) {
        @uniprots = sort @uniprots;
        @uniprots = uniq @uniprots;
        foreach my $_uniprot (@uniprots) { 
            @pmids = sort @pmids;
            @pmids = uniq @pmids;
            my $_pmid_text = join "; ", @pmids;
            print "$_uniprot\t$_pmid_text\t\n";
        }
    }
    @uniprots = ();
    @pmids    = ();
    return;
}

