#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $list       = q{};
my $gff        = q{};
my $map_target = q{};

$list       = $ARGV[0] if $ARGV[0];
$gff        = $ARGV[1] if $ARGV[1];
$map_target = $ARGV[2] if $ARGV[2];

my @genes = ();
my $data_ref;

my $header = "Gene\tMap_genome\tMap_status";

if ( (! $list ) or (! $gff ) or (! $map_target ) ) {
    die "Format: lifton_gff_genes_12feb2026.pl [input gene list] [GFF3] [mapped genome] > [table of mapped genes]\n";
}

open my $LIST, '<', $list;
while ( my $gene = <$LIST> ) {
    chomp $gene;
    if ( $gene =~ /\A \S+ \z/xms ) {
        $data_ref->{'gene'}->{$gene}->{'observed'} = 1;
        push @genes, $gene;
    }
    else {
        die "From gene list $list, cannot parse gene name: $gene\n";
    }
}
close $LIST;

open my $GFF, '<', $gff;
while ( my $input = <$GFF> ) {
    chomp $input;
    if ( $input =~ /\A [^\t]* \t [^\t]* \t gene \t /xms ) {
         if ( $input =~ /\A [^\t]* \t [^\t]* \t gene \t \d+ \t \d+ \t [^\t]* \t [+|-] \t [^\t]* \t ID[=](\S+?)[;]extra_copy_number[=]\d+; /xms ) {
             my $geneid = $1;
             my $gene   = q{};
             
             if ( $geneid !~ /\A \S+_\d+ \z/xms ) {
                 die "From GFF $gff, nonsuffixed yet multimapped gene name in: $input\n";
             }

             if ( $geneid =~ /\A (\S+) _\d+ \z/xms ) {
                 $gene = $1;
             }

             if ( exists $data_ref->{'gene'}->{$gene}->{'observed'} ) {
                 $data_ref->{'gene'}->{$gene}->{'status'} = "Multimapped";
             }
             else {
                 die "From GFF $gff, gene not in the original gene list $list: $gene\n";
             }
         }
         elsif ( $input =~ /\A [^\t]* \t [^\t]* \t gene \t \d+ \t \d+ \t [^\t]* \t [+|-] \t [^\t]* \t ID[=](\S+?)[;] /xms ) {
             my $gene = $1;

             if ( $gene =~ /_\d+\z/xms ) {
                 die "From GFF $gff, suffixed yet uniquely mapped gene name in: $input\n";
             }

             if ( exists $data_ref->{'gene'}->{$gene}->{'observed'} ) {
                 if (! exists $data_ref->{'gene'}->{$gene}->{'status'} ) {
                     $data_ref->{'gene'}->{$gene}->{'status'} = "Unique";
                 }
             }
             else {
                 die "From GFF $gff, gene not in the original gene list $list: $gene\n";
             }
         }
         else {
             die "From GFF $gff, cannot parse gene name in: $input\n";
         }
    }
}
close $GFF;

@genes = uniq(@genes);

foreach my $gene (@genes) {
    my $status = "Unmapped";

    if ( exists $data_ref->{'gene'}->{$gene}->{'status'} ) {
        $status = $data_ref->{'gene'}->{$gene}->{'status'};
    }

    print "$header\n" if $header;
    $header = q{};
    print "$gene\t$map_target\t$status\n";
}
