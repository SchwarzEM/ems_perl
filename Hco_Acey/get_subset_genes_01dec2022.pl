#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $acceptable = q{};
my $raw_genes  = q{};
my %seen       = ();

$acceptable  = $ARGV[0] if $ARGV[0];
$raw_genes = $ARGV[1] if $ARGV[1];

if ( (! -e $acceptable) or (! -e $raw_genes) ) {
    die "Format: get_subset_genes_01dec2022.pl",
        " [list of acceptable genes] [unfiltered list of genes]",
        " > [filtered list of genes]",
        "\n";
}

open my $ACCEPTABLE, '<', $acceptable;
while (my $input = <$ACCEPTABLE>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \z/xms ) {
        my $ok_gene = $1;
        $seen{$ok_gene} = 1;
    }
    else {
        die "From list of acceptable genes $acceptable, cannot parse: $input\n";
    }
}
close $ACCEPTABLE;

open my $RAW_GENES, '<', $raw_genes;
while (my $input = <$RAW_GENES>) {
    chomp $input; 
    if ( $input =~ /\A (\S+) \z/xms ) {
        if ( exists $seen{$input} ) {
            print "$input\n";
        }
    }
    else {
        die "From raw_genes list $raw_genes, cannot parse: $input\n";
    }
}
close $RAW_GENES;

