#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $infile    = q{};
my $key_taxon = q{};
my $header    = "Gene\tPangene\tPg_PASS_taxa\tPg_array_taxa\tPg_NSOrtho_taxa";
my @pangenes  = ();

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
        my @genes = sort keys %{ $data_ref->{'pangene'}->{$pangene}->{'taxon'}->{'Aroian'}->{'flag'}->{'PASS'}->{'gene'} };
        foreach my $gene (@genes) {
            my @pass_taxa = sort keys %{ $data_ref->{'pangene'}->{$pangene}->{'flag'}->{'PASS'}->{'taxon'} };
            my $pass_taxon_text = join '; ', @pass_taxa;

            my @array_taxa = sort keys %{ $data_ref->{'pangene'}->{$pangene}->{'flag'}->{'array'}->{'taxon'} };
            my $array_taxon_text = join '; ', @array_taxa;

            my @nsortho_taxa = sort keys %{ $data_ref->{'pangene'}->{$pangene}->{'flag'}->{'NSOrtho'}->{'taxon'} };
            my $nsortho_taxon_text = join '; ', @nsortho_taxa;

            if ( $header ) {
                print "$header\n";
                $header = q{};
            }

            print "$gene\t$pangene\t$pass_taxon_text\t$array_taxon_text\t$nsortho_taxon_text\n";
        }
    }
}

