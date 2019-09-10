#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Getopt::Long;

my $ofind  = q{};
my @taxa   = ();
my @counts = ();
my $help;

my $data_ref;

GetOptions ( 'ofind=s'     => \$ofind,
             'taxa=s{,}'   => \@taxa, 
             'counts=s{,}' => \@counts,
             'help'        => \$help, );

if ( $help or (! $ofind) or (! @taxa ) or (! @counts ) ) {
    die "Format: add_aa_sum_columns_09sep2019.pl\n",
        "    --ofind|-o   [gene table containing one OrthoFinder entry per line]\n",
        "    --taxa|-t    [taxa for which to provide amino acid sums per OrthoFinder entry]\n",
        "    --counts|-c  [list of gene-to-max_aa tables; order must match taxa]\n",
        "    --help|-h    [print this message]\n",
        ;
}

my $j = @taxa;
my $k = @counts;
if ( $j != $k ) {
    die "Unequal counts in taxa versus counts\n";
}
$j--;

foreach my $i (0..$j) {
    my $taxon     = $taxa[$i];
    my $countfile = $counts[$i];
    open my $COUNTFILE, '<', $countfile;
    while (my $input = <$COUNTFILE>) {
        chomp $input;
        if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
            my $gene     = $1;
            my $aa_count = $2;
            if ( ( $gene ne 'Gene' ) and ( $aa_count =~ /\d+/xms ) ) {
                if ( exists $data_ref->{'taxon'}->{$taxon}->{'gene'}->{$gene} ) {
                    die "Redundant values for aa count of taxon $taxon, gene $gene: ",
                        "$data_ref->{'taxon'}->{$taxon}->{'gene'}->{$gene} and $aa_count\n";
                }
                $data_ref->{'taxon'}->{$taxon}->{'gene'}->{$gene} = $aa_count;
            }
            elsif ( $gene ne 'Gene' ) {
                 die "In count file $countfile, cannot decipher: $input\n";
            }
        }
        else {
            die "In count file $countfile, cannot parse: $input\n";
        }
    }
    close $COUNTFILE;
}

open my $OFIND, '<', $ofind;
while (my $input = <$OFIND>) {
    chomp $input;
    my $output = $input;
    if ( $input =~ /\A Gene \t (?:OrthoMCL|OrthoFinder) \z/xms ) {
        foreach my $taxon (@taxa) {
            my $header = "$taxon.count";
            $output    = "$output\t$header";
        }
        print "$output\n";
    }
    elsif ( $input =~ /\A (\S+) \t \S+\(\d+\sgenes,\d+\staxa\): \s+ (\S.+\S) \z/xms ) { 
        my $gene       = $1;
        my $ogene_text = $2;
        my @ogenes = split /\s+/, $ogene_text;

        # zero this out
        delete $data_ref->{'ofind'};

        # get aa numbers, add them up for each taxon of the ofind group
        foreach my $ogene (@ogenes) {
            if ( $ogene =~ /(\S+) \( (\S+) \)/xms ) {
                my $og_id  = $1;
                my $og_tax = $2;
                if (! exists $data_ref->{'taxon'}->{$og_tax}->{'gene'}->{$og_id} ) {
                    die "Cannot get aa count for gene $ogene in: $input\n";
                }
                my $og_aa = $data_ref->{'taxon'}->{$og_tax}->{'gene'}->{$og_id};
                $data_ref->{'ofind'}->{'taxon'}->{$og_tax} += $og_aa;
            }
            else {
                die "In orthology file $ofind, cannot parse gene $ogene in: $input\n";
            }
        }

        foreach my $taxon (@taxa) {
            my $aa_total = $data_ref->{'ofind'}->{'taxon'}->{$taxon};
            $output     = "$output\t$aa_total";
        }
        print "$output\n";
    }
    else {
        die "In orthology file $ofind, cannot parse: $input\n";
    }
}
close $OFIND;
