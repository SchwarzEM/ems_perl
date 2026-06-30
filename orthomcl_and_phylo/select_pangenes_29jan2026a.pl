#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $infile      = q{};
my $key_taxon   = q{};
my $header      = "Gene\tPangenes\tPg_PASS_taxa\tPg_array_taxa\tPg_NSOrtho_taxa";
my @pangenes    = ();
my @final_genes = ();

my $data_ref;

$infile    = $ARGV[0] if $ARGV[0];
$key_taxon = $ARGV[1] if $ARGV[1];

if ( (! $infile ) or (! $key_taxon ) ) {
    die "Format: select_pangenes_23jan2026a.pl [infile] [key taxon] > [key taxon syntelogs]\n";
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
        $data_ref->{'pangene'}->{$pangene}->{'taxon'}->{$taxon}->{'flag'}->{$flag}->{'gene'}->{$gene} = 1;
        $data_ref->{'pangene'}->{$pangene}->{'flag'}->{$flag}->{'taxon'}->{$taxon}->{'gene'}->{$gene} = 1;
        push @pangenes, $pangene;
    }
}
close $INFILE;

@pangenes = uniq(@pangenes);

foreach my $pangene (@pangenes) {
    if ( exists $data_ref->{'pangene'}->{$pangene}->{'taxon'}->{$key_taxon}->{'flag'}->{'PASS'}->{'gene'} ) {
        my @genes = sort keys %{ $data_ref->{'pangene'}->{$pangene}->{'taxon'}->{$key_taxon}->{'flag'}->{'PASS'}->{'gene'} };
        foreach my $gene (@genes) {
            $data_ref->{'gene'}->{$gene}->{'pangene'}->{$pangene} = 1;

            my @pass_taxa = sort keys %{ $data_ref->{'pangene'}->{$pangene}->{'flag'}->{'PASS'}->{'taxon'} };
            foreach my $pass_taxon (@pass_taxa) {
                $data_ref->{'gene'}->{$gene}->{'pass_taxon'}->{$pass_taxon} = 1;
            }

            my @array_taxa = sort keys %{ $data_ref->{'pangene'}->{$pangene}->{'flag'}->{'array'}->{'taxon'} };
            foreach my $array_taxon (@array_taxa) {
                $data_ref->{'gene'}->{$gene}->{'array_taxon'}->{$array_taxon} = 1;
            }

            my @nsortho_taxa = sort keys %{ $data_ref->{'pangene'}->{$pangene}->{'flag'}->{'NSOrtho'}->{'taxon'} };
            foreach my $nsortho_taxon (@nsortho_taxa) {
                $data_ref->{'gene'}->{$gene}->{'nsortho_taxon'}->{$nsortho_taxon} = 1;
            }
            push @final_genes, $gene;
        }
    }
}

@final_genes = uniq(@final_genes);

foreach my $gene (@final_genes) { 
    my @g_pangenes   = sort keys %{ $data_ref->{'gene'}->{$gene}->{'pangene'} };
    my @pass_taxa    = sort keys %{ $data_ref->{'gene'}->{$gene}->{'pass_taxon'} };
    my @array_taxa   = sort keys %{ $data_ref->{'gene'}->{$gene}->{'array_taxon'} };
    my @nsortho_taxa = sort keys %{ $data_ref->{'gene'}->{$gene}->{'nsortho_taxon'} };

    my $pangene_text       = join '; ', @g_pangenes;
    my $pass_taxon_text    = join '; ', @pass_taxa;
    my $array_taxon_text   = join '; ', @array_taxa;
    my $nsortho_taxon_text = join '; ', @nsortho_taxa;

    if ( $header ) {
        print "$header\n";
        $header = q{};
    }

    print "$gene\t$pangene_text\t$pass_taxon_text\t$array_taxon_text\t$nsortho_taxon_text\n";
}
