#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $gene_list = q{};
my $ofinds    = q{};
my $blastps   = q{};
my $orphans   = q{};

$gene_list = $ARGV[0] if $ARGV[0];
$ofinds    = $ARGV[1] if $ARGV[1];
$blastps   = $ARGV[2] if $ARGV[2];
$orphans   = $ARGV[3] if $ARGV[3];

my @gene_list = ();
my $data_ref;

my $header = "Gene\tOrphan\/non-orphan";

open my $LIST, '<', $gene_list;
while (my $input = <$LIST>) {
    chomp $input;
    if ( $input ne 'Gene' ) {
        if ( $input =~ /\A \S+ \z/xms ) { 
            push @gene_list, $input;
            $data_ref->{'listed_gene'}->{$input} = 1;
        }
        else { 
            die "In gene list $gene_list, cannot parse gene name \"$input\"\n";
        }
    }
}
close $LIST;

@gene_list = uniq(@gene_list);

open my $OFINDS, '<', $ofinds;
while (my $input = <$OFINDS>) {
    chomp $input;
    if ( $input =~ / \A (\S+) \t OG\d+ \( \d+ [ ] genes , \d+ [ ] taxa \) [:] [ ] (.+) \z/xms ) { 
        my $gene          = $1;
        my $ofind_summary = $2;

        if (! exists $data_ref->{'listed_gene'}->{$gene} ) {
            die "In OrthoFinder summary $ofinds, cannot recognize gene \"$gene\" as listed gene: $input\n";
        }

        my @ofind_taxa = split '; ', $ofind_summary;
        foreach my $ofind_taxon_text (@ofind_taxa) {
            if ( $ofind_taxon_text =~ /\A (\S+) [ ] \( \d+ [ ] g\. \) \z/xms ) {
                my $ofind_taxon = $1;
                if ( $ofind_taxon ne 'nigoni' ) {
                    $data_ref->{'gene'}->{$gene}->{'ofind_taxon'}->{$ofind_taxon} = 1;
                }
            }
            else { 
                die "In OrthoFinder summary $ofinds, cannot parse taxa \"$ofind_taxon_text\" in: $input\n";
            }
        }
    }
    elsif (  ( $input !~ /\A Gene \t /xms ) and ( $input !~ /\A \S+ \t \z/xms ) ) { 
        die "In OrthoFinder summary $ofinds, cannot parse text line: $input\n";
    }
}
close $OFINDS;

open my $BLASTPS, '<', $blastps;
while (my $input = <$BLASTPS>) {
    chomp $input;
    if ( $input =~ / \A (\S+) \t ([^\t]+) \z/xms ) { 
        my $gene         = $1;
        my $blastps_text = $2;

        if (! exists $data_ref->{'listed_gene'}->{$gene} ) {
            die "In BlastP summary $blastps, cannot recognize gene \"$gene\" as listed gene: $input\n";
        }

        my @blastps_taxa = split '; ', $blastps_text;
        foreach my $blastps_taxon (@blastps_taxa) {
            if ( $blastps_taxon ne 'nigoni' ) {
                $data_ref->{'gene'}->{$gene}->{'blastp_taxon'}->{$blastps_taxon} = 1;
            }
        }
    }
    else { 
        die "From BlastP summary $blastps, cannot parse: $input\n";
    }
}
close $BLASTPS;

open my $ORPHANS, '<', $orphans;
while (my $input = <$ORPHANS>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \z/xms ) { 
        my $gene = $1;

        if (! exists $data_ref->{'listed_gene'}->{$gene} ) {
            die "In orphan gene list $orphans, cannot recognize as listed gene: $input\n";
        }
        if ( exists $data_ref->{'gene'}->{$gene}->{'ofind_taxon'} ) {
            die "Gene $gene is listed as both being an orphan and as having OrthoFinder non-nigoni taxa!\n";
        }
        $data_ref->{'gene'}->{$gene}->{'orphan_gene'} = 1;
    }
}
close $ORPHANS;

foreach my $listed_gene (@gene_list) {
    # Do this exactly once, at the start:
    print "$header\n" if $header;
    $header = q{};

    if ( exists $data_ref->{'gene'}->{$listed_gene}->{'blastp_taxon'} ) { 
        my @blastp_taxa      = sort keys %{ $data_ref->{'gene'}->{$listed_gene}->{'blastp_taxon'} };
        my $blastp_taxa_text = join '; ', @blastp_taxa;
        print "$listed_gene\tBlastP [$blastp_taxa_text]\n";
    }
    elsif ( exists $data_ref->{'gene'}->{$listed_gene}->{'orphan_gene'} ) { 
        print "$listed_gene\tOrphan\n";
    }
    elsif ( exists $data_ref->{'gene'}->{$listed_gene}->{'ofind_taxon'} ) { 
        my @ofind_taxa = sort keys %{ $data_ref->{'gene'}->{$listed_gene}->{'ofind_taxon'} };

        my $ofind_taxa_text = join '; ', @ofind_taxa;
        print "$listed_gene\tOrthoFinder [$ofind_taxa_text]\n";   
    }
    else { 
        die "Cannot determine homology status of gene: $listed_gene\n";
    }
}

