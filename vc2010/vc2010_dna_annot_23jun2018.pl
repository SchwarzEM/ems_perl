#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $one_nt_genes   = $ARGV[0] if $ARGV[0];
my $one_exon_genes = $ARGV[1] if $ARGV[1];
my $all_exon_genes = $ARGV[2] if $ARGV[2];

my $header = "Gene\tVC2010-spec.";

my $data_ref;

open my $NT_1, '<', $one_nt_genes;
while (my $input = <$NT_1>) {
    chomp $input;
    if ( $input =~ /\A \S+ \z/xms ) {
        $data_ref->{'gene'}->{$input} = '1+_nt';
    }
}
close $NT_1;

open my $EXON, '<', $one_exon_genes;
while (my $input = <$EXON>) {
    chomp $input; 
    if ( $input	=~ /\A \S+ \z/xms ) {
       	$data_ref->{'gene'}->{$input} =	'1+_exons';
    }
}
close $EXON;

open my $ALL, '<', $all_exon_genes;
while (my $input = <$ALL>) {
    chomp $input; 
    if ( $input	=~ /\A \S+ \z/xms ) {
       	$data_ref->{'gene'}->{$input} =	'All_exons';
    }
}
close $ALL;

my @genes = sort keys %{ $data_ref->{'gene'} };

foreach my $gene (@genes) {
    print "$header\n" if $header;
    $header = q{};

    my $annot = $data_ref->{'gene'}->{$gene};
    print "$gene\t$annot\n";
}
