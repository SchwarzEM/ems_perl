#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $map   = q{};
my $annot = q{};

my $data_ref;

$map   = $ARGV[0] if $ARGV[0];
$annot = $ARGV[1] if $ARGV[1];

if ( (! $map ) or (! $annot ) ) {
    die "Format: washu_v1_to_acey_v2.1_03jul2023.pl [mapping of WashU to Acey genes] [annotation table] > [remapped annotation] ;\n";
}

open my $MAP, '<', $map;
while (my $input = <$MAP>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S.*) \z/xms ) {
        my $acey = $1;
        my $data = $2;
        if ( $acey ne 'Gene' ) {
            my @washus = split /; /, $data;
            foreach my $washu (@washus) {
                if ( $washu =~ /\S/xms ) {
                    $data_ref->{'washu'}->{$washu}->{'acey'}->{$acey} = 1;
                }
            }
        }
    }
    else {
        die "From map file $map, cannot parse: $input\n";
    }
}
close $MAP;

open my $ANNOT, '<', $annot;
while (my $input = <$ANNOT>) {

    chomp $input;
    if ( $input =~ /\A (\S+) \t (.+) \z/xms ) {
        my $gene  = $1;
        my $annot = $2;
        if ( exists $data_ref->{'washu'}->{$gene} ) {
            my @mapped_genes = sort keys %{ $data_ref->{'washu'}->{$gene}->{'acey'} };
            foreach my $mapped_gene (@mapped_genes) {
                print "$mapped_gene\t$annot\n";
            }
        }
        else {
            print "$gene\t$annot\n";
        }
    }
    else {   
        die "From annot file $annot, cannot parse: $input\n";
    }
}
close $ANNOT;

