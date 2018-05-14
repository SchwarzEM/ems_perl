#!/usr/bin/env perl

# link_wbgenes2protmots.pl -- Erich Schwarz <emsch@its.caltech.edu>, 4/18/2008.
# Purpose: given gene-ID and prot-motif lists, summarize gene-centrically.

use strict;
use warnings;

my %wbgeneids_ref  = ();
my %wbgenemots_ref = ();

while ( my $input = <> ) {
    chomp $input;

    # Decide which file is being read solely by line format (1 or 2 below):

    # Examples of format 1:
    # WBGene00000262  F54B11.6        bra-1
    # WBGene00000263  F23H11.5

    if (
        $input =~ / \A 
                     ( WBGene\d+ ) 
                     \t 
                     # Ensure that CGC-less gets captured:
                     ( \S+ \t .* ) 
                     \z
                   /xms
      )
    {
        my $wbgene     = $1;
        my $wbgene_ids = $2;
        if ( $wbgene_ids =~ / \A (\S+) /xms ) {
            my $cds = $1;
            if ( $cds !~ /\./xms ) {
                die "$cds in $input is misformatted CDS name.\n";
            }
            $wbgeneids_ref{$wbgene}->{'cds'} = $cds;
        }
        if ( $wbgene_ids =~ / \A \S+ \t (\S+) \z /xms ) {
            my $cgc = $1;
            if ( $cgc !~ /-/xms ) {
                die "$cgc in $input is misformatted CGC name.\n";
            }
            $wbgeneids_ref{$wbgene}->{'cgc'} = $cgc;
        }
    }

    # Examples of format 2:
    # WP:CE05930      InP_met_003784          WBGene00000262
    # WP:CE05930      KOG3612 PHD Zn-finger protein   WBGene00000262
    # WP:CE05930      OMpre_WH000540          WBGene00000262

    if (
        $input =~ / \A
                     WP:CE\d+ 
                     \t
                     ( \S+ ) 
                     \t 
                     ( [^\t]+ )
                     \t
                     ( WBGene\d+ )
                     \z
                  /xms
      )
    {
        my $kog_or_mot  = $1;
        my $description = $2;
        my $wbgene      = $3;
        $kog_or_mot = "$kog_or_mot [$description]";

        # Prevent redundant descriptors:
        $wbgenemots_ref{$wbgene}->{$kog_or_mot} = 1;
    }
}

foreach my $wbgene ( sort keys %wbgeneids_ref ) {

    # Print a sensible gene ID column:
    print $wbgene;
    print q{|};
    print $wbgeneids_ref{$wbgene}->{'cds'};
    if ( $wbgeneids_ref{$wbgene}->{'cgc'} ) {
        print q{|};
        print $wbgeneids_ref{$wbgene}->{'cgc'};
    }
    print "\t";

    # Print descriptions, if they exist for the gene:
    my @descriptors;
    if ( exists $wbgenemots_ref{$wbgene} ) {
        @descriptors = sort keys %{ $wbgenemots_ref{$wbgene} };
    }
    my $descriptor;
    if (@descriptors) {
        $descriptor = join '; ', @descriptors;
        print $descriptor;
    }

    # End the line.
    print "\n";
}

