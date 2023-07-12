#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $prot_table = q{};
$prot_table    = $ARGV[0] if $ARGV[0];

my @seqs = ();

if (! $prot_table ) {
    die "Format: non_redund_aug_txs_11jul2023.pl [table from prots_to_seqids_11jul2023.pl] > [list of nonredundant txs within genes]\n";
}

open my $PROT_TABLE, '<', $prot_table;
while ( my $input = <$PROT_TABLE> ) {
    chomp $input;
    if ( $input =~ /\A \S+ \z/xms ) {
        push @seqs, $input;
    }
    elsif ( $input =~ /[;][ ]/xms ) {
        my @txs = split /; /, $input;
        my $data_ref;
        foreach my $tx (@txs) {
            # Sample input: Nippo_chrV.g13264.t1 
            if ( $tx =~ /\A (\S+) \. t (\d+) \z/xms ) {
                my $gene  = $1;
                my $tx_no = $2;
                $data_ref->{'gene'}->{$gene}->{'tx_no'}->{$tx_no}->{'tx'} = $tx;
            }
            else {
                die "From prot_table $prot_table, cannot parse tx \"$tx\" in: $input\n";
            }
        }
        my @genes = sort keys %{ $data_ref->{'gene'} };
        foreach my $gene (@genes) {
            if ( exists $data_ref->{'gene'}->{$gene}->{'tx_no'} ) {
                my @tx_nos = sort { $a <=> $b } keys %{ $data_ref->{'gene'}->{$gene}->{'tx_no'} };
                my $tx_no  = $tx_nos[0];
                my $tx_rep = q{};
                if ( exists $data_ref->{'gene'}->{$gene}->{'tx_no'}->{$tx_no}->{'tx'} ) {
                    $tx_rep = $data_ref->{'gene'}->{$gene}->{'tx_no'}->{$tx_no}->{'tx'};
                }
                else {
                    die "Cannot map gene $gene and tx_no $tx_no to transcript in: $input\n";
                }
                push @seqs, $tx_rep;
            }
            else {
                die "Failed to map gene $gene to transcript numbers\n";
            }
        }
    }
}
close $PROT_TABLE;

foreach my $seq (@seqs) {
    print "$seq\n";
}
