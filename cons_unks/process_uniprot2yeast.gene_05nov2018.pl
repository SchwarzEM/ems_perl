#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $infile = q{};
my $prefix = q{};

$infile = $ARGV[0] if $ARGV[0];
$prefix = $ARGV[1] if $ARGV[1];

open my $INFILE, '<', $infile;
while (my $input = <$INFILE>) {
    chomp $input;
    if ( $input =~ /\A [a-z]+\|(\S+)\| [^\t]* \t (.*) \z/xms ) {
        my $uniprot        = $1;
        my $raw_gene_data  = $2;
        my @raw_gene_texts = split /\t/, $raw_gene_data;
        my @gene_names     = ();

        foreach my $raw_gene_text (@raw_gene_texts) {
            if ( $raw_gene_text =~ /\A \S+ \s+ (\S+) .+ SGDID:(S\d+) /xms ) {
                my $gene_name = $1;
                my $gene_id   = $2;
                my $gene      = "$prefix|$gene_id|$gene_name";
                push @gene_names, $gene;
            }
            else {
                die "Can't parse gene info in: $raw_gene_text\n";
            }
        }
        @gene_names = sort @gene_names;
        @gene_names = uniq @gene_names;
        my $gene_list = join '; ', @gene_names;
        print "$uniprot\t$gene_list\n";
    }
}
close $INFILE;

