#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my $infile = q{};
$infile    = $ARGV[0] if $ARGV[0];

my $header = "Gene\tMax_diff_ab\tAssoc_p-value\tAssoc_diff_ab\tMin_p-value\n";

if (! $infile ) {
     die "Format: map_ichip_to_acey_genes_22may2021.pl [data file] > [gene data table]\n";
}

open my $INFILE, '<', $infile;
while (my $input = <$INFILE>) {
    if ( $input =~ /\A \S+ \t [^\t]* \t ([^\t]+) (?: \t [^\t]*){34} \t ([^\t]+) \t ([^\t]+) \t /xms ) {
        my $acey_text = $1;
        my $diff_ab   = $2;
        my $p_value   = $3;
        my @acey_homs = split /[;] /, $acey_text;
        foreach my $acey_hom (@acey_homs) {
            if ( $acey_hom =~ /\A (\S+) \(a_ceylanicum\) \z/xms ) {
                my $acey_gene = $1;
                $data_ref->{'Acey_gene'}->{$acey_gene}->{'diff_ab'}->{$diff_ab}->{'assoc_p_value'}->{$p_value} = 1;
                $data_ref->{'Acey_gene'}->{$acey_gene}->{'p_value'}->{$p_value}->{'assoc_diff_ab'}->{$diff_ab} = 1;
            }
            elsif ( $acey_hom ne 'Acey_homologs' ) {
                die "Cannot parse Acey gene: $acey_hom\n";
            }
        }
    }
}

my @acey_genes = sort keys %{ $data_ref->{'Acey_gene'} };
foreach my $acey_gene (@acey_genes) {
    my @diff_abs = sort { $b <=> $a } keys %{ $data_ref->{'Acey_gene'}->{$acey_gene}->{'diff_ab'} };
    my @p_values = sort { $a <=> $b } keys %{ $data_ref->{'Acey_gene'}->{$acey_gene}->{'p_value'} };

    my $max_diff_ab    = $diff_abs[0];
    my @assoc_p_values = sort { $a <=> $b } keys %{ $data_ref->{'Acey_gene'}->{$acey_gene}->{'diff_ab'}->{$max_diff_ab}->{'assoc_p_value'} };
    my $assoc_p_value  = $assoc_p_values[0];

    my $min_p_value    = $p_values[0];
    my @assoc_diff_abs = sort { $b <=> $a } keys %{ $data_ref->{'Acey_gene'}->{$acey_gene}->{'p_value'}->{$min_p_value}->{'assoc_diff_ab'} };
    my $assoc_diff_ab  = $assoc_diff_abs[0];
    
    print "$header\n" if $header;
    $header = q{};

    print "$acey_gene\t$max_diff_ab\t$assoc_p_value\t$assoc_diff_ab\t$min_p_value\n";
}

