#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $prefix = q{};
my $infile = q{};

my $header     = "Gene\tclust_group";
my %gene2clust = ();

$prefix = $ARGV[0] if $ARGV[0];
$infile = $ARGV[1] if $ARGV[1];

foreach my $i (0..13) {
    my $j = ($i + 1);
    my @inputs = `cat $infile | cut -f $j ;`;
    my $clust_name = q{};
    foreach my $input (@inputs) {
        chomp $input;

        # Sample input, first line:  C0 (5753 genes)
        if ( $input =~ /\A C\d+ \s+ \( \d+ \s+ genes \) \s* \z/xms) {
            $clust_name = $input;

            $clust_name =~ s/\s*\z//;
            $clust_name =~ s/\s+/_/g;
            $clust_name =~ s/\(/_/g;
            $clust_name =~ s/\)/_/g;
            $clust_name =~ s/[_]+/_/g;
            $clust_name =~ s/_\z//;
            $clust_name = "$prefix.$clust_name";
        }

        # Don't record the 'Genes' subheader; also don't record empty lines!
        elsif ( ( $input ne 'Genes' ) and ( $input =~ /\S/xms ) ) {
            $gene2clust{$input} = $clust_name;
        }
    }
}

my @genes = sort keys %gene2clust;
foreach my $gene (@genes) {
    # Do this only once, at the start:
    print "$header\n" if $header;
    $header = q{};

    my $clust_name = $gene2clust{$gene};
    print "$gene\t$clust_name\n";
}


