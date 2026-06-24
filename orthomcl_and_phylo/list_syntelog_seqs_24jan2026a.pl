#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $infile    = q{};
my $key_taxon = q{};
my @pangenes  = ();

my $data_ref;

$infile    = $ARGV[0] if $ARGV[0];
$key_taxon = $ARGV[1] if $ARGV[1];

if ( (! $infile ) or (! $key_taxon ) ) {
    die "Format: select_pangenes_23jan2026a.pl [infile] [key taxon] > [pangene name; taxon|syntelog list]\n";
}

open my $INFILE, '<', $infile;
while ( my $input = <$INFILE> ) {
    chomp $input;

    # Sample input:
    # ofID    pgID    interpChr       interpOrd       pgRepID genome  og      flag    id      chr     start   end     ord
    # 0_368   1       Necator_chrI    1       0_368   Aroian  35669   PASS    Necator_chrI.g137       Necator_chrI    13892   14485   1

    if ( ( $input !~ /\A ofID \t/xms ) and 
         ( $input =~ /\A \S+ \t (\S+) \t (?: \S+ \t){3} (\S+) \t (?: \S+ \t) (\S+) \t (\S+) \t .* \z/xms ) ) {
        my $pangene = $1;
        my $taxon   = $2;
        my $flag    = $3;
        my $gene    = $4;
        $data_ref->{'pangene'}->{$pangene}->{'flag'}->{$flag}->{'taxon'}->{$taxon}->{'gene'}->{$gene} = 1;
        push @pangenes, $pangene;
    }
}
close $INFILE;

@pangenes = uniq(@pangenes);

foreach my $pangene (@pangenes) {
    if ( exists $data_ref->{'pangene'}->{$pangene}->{'flag'}->{'PASS'}->{'taxon'}->{$key_taxon}->{'gene'} ) {
        my @taxa = sort keys %{ $data_ref->{'pangene'}->{$pangene}->{'flag'}->{'PASS'}->{'taxon'} };
        my @tax_genes = ();
        foreach my $taxon (@taxa) {
            my @genes = sort keys %{ $data_ref->{'pangene'}->{$pangene}->{'flag'}->{'PASS'}->{'taxon'}->{$taxon}->{'gene'} };
            foreach my $gene (@genes) {
                my $tax_gene = "$taxon|$gene";
                push @tax_genes, $tax_gene;
            }
        }
        @tax_genes = uniq(@tax_genes);
        my $tax_gene_text = join '; ', @tax_genes;
        print "$pangene \t$tax_gene_text\n";
    }
}

