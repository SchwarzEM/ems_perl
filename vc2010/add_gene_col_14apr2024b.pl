#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $cds2g_tab = q{};
my $cds_tab   = q{};

$cds2g_tab = $ARGV[0] if $ARGV[0];
$cds_tab   = $ARGV[1] if $ARGV[1];

my %cds2gene = ();

if ( (! $cds2g_tab ) or (! $cds_tab ) ) {
    die "Format: add_gene_col_14apr2024.pl [cds2gene table] [cds + something table] > [cds + gene + something table]\n";
}

open my $CDS2G_TAB, '<', $cds2g_tab;
while ( my $input = <$CDS2G_TAB> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $cds  = $1;
        my $gene = $2;
        if ( exists $cds2gene{$cds} ) {
            die "Redundant mapping of CDS $cds to both $cds2gene{$cds} and $gene\n";
        }
        $cds2gene{$cds} = $gene;
    }
    else {
        die "From cgc2gene table $cds2g_tab, cannot parse: $input\n";
    }
}
close $CDS2G_TAB;

open my $CDS_TAB, '<', $cds_tab;
while ( my $input = <$CDS_TAB> ) {
    chomp $input; 
    if ( $input =~ /\A (\S+) \t (.+) \z/xms ) { 
        my $cds   = $1;
        my $annot = $2;
        my $gene  = q{};
        if (! exists $cds2gene{$cds} ) {
            my $alt_cds = $cds;
            $alt_cds =~ s/_\d+\z//;

            if (! exists $cds2gene{$alt_cds} ) {
                die "Cannot map CDS $cds to a gene\n";
            }
            $gene = $cds2gene{$alt_cds};

            my $suffix = q{};
            if ( $cds =~ /\A \S+ (_\d+) \z/xms ) {
                $suffix = $1;
            }
            else {
                die "Cannot extract suffix from CDS $cds\n";
            }

            # ncRNA sequence names do not need to be stripped of the extra '.1+' terms for gffread.
            # The gene name still does need to be revised human-readably.
            $gene = $gene . '.extra_copy' . $suffix;
        }
        else {
            $gene = $cds2gene{$cds};
        }
        print "$cds\t$gene\t$annot\n";
    }
    else {
        die "From CDS-annotation table $cds_tab, cannot parse: $input\n";
    }
}
close $CDS_TAB;

