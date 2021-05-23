#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my $infile = q{};
$infile    = $ARGV[0] if $ARGV[0];

my $header = "Gene" 
             . "\t"

             . "Max_diff_all" 
             . "\t" 
             . "Assoc_p-value_all" 
             . "\t" 

             . "Min_p-value_all"
             . "\t"
             . "Assoc_diff_all" 
             . "\t" 

             . "Max_diff_hw"
             . "\t"
             . "Assoc_p-value_hw"
             . "\t"

             . "Min_p-value_hw"
             . "\t"
             . "Assoc_diff_hw"
             ;

if (! $infile ) {
     die "Format: map_ichip_to_acey_genes_22may2021.pl [data file] > [gene data table]\n";
}

open my $INFILE, '<', $infile;
while (my $input = <$INFILE>) {
    if ( $input =~ /\A \S+ \t [^\t]* \t ([^\t]+) (?: \t [^\t]*){34} \t ([^\t]+) \t ([^\t]+) (?: \t [^\t]*){3} \t ([^\t]+) \t ([^\t]+) /xms ) {
        my $acey_text   = $1;
        my $diff_all    = $2;
        my $p_value_all = $3;
        my $diff_hw     = $4;
        my $p_value_hw  = $5;
        my @acey_homs = split /[;] /, $acey_text;
        foreach my $acey_hom (@acey_homs) {
            if ( $acey_hom =~ /\A (\S+) \(a_ceylanicum\) \z/xms ) {
                my $acey_gene = $1;
                $data_ref->{'Acey_gene'}->{$acey_gene}->{'diff_all'}->{$diff_all}->{'assoc_p_value_all'}->{$p_value_all} = 1;
                $data_ref->{'Acey_gene'}->{$acey_gene}->{'p_value_all'}->{$p_value_all}->{'assoc_diff_all'}->{$diff_all} = 1;
                $data_ref->{'Acey_gene'}->{$acey_gene}->{'diff_hw'}->{$diff_hw}->{'assoc_p_value_hw'}->{$p_value_hw} = 1;
                $data_ref->{'Acey_gene'}->{$acey_gene}->{'p_value_hw'}->{$p_value_hw}->{'assoc_diff_hw'}->{$diff_hw} = 1;
            }
            elsif ( $acey_hom ne 'Acey_homologs' ) {
                die "Cannot parse Acey gene: $acey_hom\n";
            }
        }
    }
}

my @acey_genes = sort keys %{ $data_ref->{'Acey_gene'} };
foreach my $acey_gene (@acey_genes) {
    my @diff_alls = sort { $b <=> $a } keys %{ $data_ref->{'Acey_gene'}->{$acey_gene}->{'diff_all'} };
    my @p_value_alls = sort { $a <=> $b } keys %{ $data_ref->{'Acey_gene'}->{$acey_gene}->{'p_value_all'} };

    my $max_diff_all    = $diff_alls[0];
    my @assoc_p_value_alls = sort { $a <=> $b } keys %{ $data_ref->{'Acey_gene'}->{$acey_gene}->{'diff_all'}->{$max_diff_all}->{'assoc_p_value_all'} };
    my $assoc_p_value_all  = $assoc_p_value_alls[0];

    my $min_p_value_all    = $p_value_alls[0];
    my @assoc_diff_alls = sort { $b <=> $a } keys %{ $data_ref->{'Acey_gene'}->{$acey_gene}->{'p_value_all'}->{$min_p_value_all}->{'assoc_diff_all'} };
    my $assoc_diff_all  = $assoc_diff_alls[0];

    my @diff_hws = sort { $b <=> $a } keys %{ $data_ref->{'Acey_gene'}->{$acey_gene}->{'diff_hw'} };
    my @p_value_hws = sort { $a <=> $b } keys %{ $data_ref->{'Acey_gene'}->{$acey_gene}->{'p_value_hw'} };

    my $max_diff_hw       = $diff_hws[0];
    my @assoc_p_value_hws = sort { $a <=> $b } keys %{ $data_ref->{'Acey_gene'}->{$acey_gene}->{'diff_hw'}->{$max_diff_hw}->{'assoc_p_value_hw'} };
    my $assoc_p_value_hw  = $assoc_p_value_hws[0];

    my $min_p_value_hw    = $p_value_hws[0];
    my @assoc_diff_hws = sort { $b <=> $a } keys %{ $data_ref->{'Acey_gene'}->{$acey_gene}->{'p_value_hw'}->{$min_p_value_hw}->{'assoc_diff_hw'} };
    my $assoc_diff_hw  = $assoc_diff_hws[0];

    print "$header\n" if $header;
    $header = q{};

    print "$acey_gene" 
          . "\t"

          . "$max_diff_all"
          . "\t"
          . "$assoc_p_value_all"
          . "\t"

          . "$min_p_value_all"
          . "\t"
          . "$assoc_diff_all"
          . "\t"

          . "$max_diff_hw"
          . "\t"
          . "$assoc_p_value_hw"
          . "\t"

          . "$min_p_value_hw"
          . "\t"
          . "$assoc_diff_hw"
          . "\n"
          ;
}

