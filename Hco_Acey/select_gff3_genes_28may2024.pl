#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $genelist = q{};
my $gff3     = q{};

my %ok_genes = ();

$genelist = $ARGV[0] if $ARGV[0];
$gff3     = $ARGV[1] if $ARGV[1];

if ( (! $genelist ) or (! $gff3 ) ) {
    die "Format: select_gff3_genes_28may2024.pl [gene list] [GFF3] > [GFF3 with listed genes]\n";
}

open my $GENELIST, '<', $genelist;
while ( my $gene = <$GENELIST> ) {
    chomp $gene;
    $ok_genes{$gene} = 1;
}
close $GENELIST;

open my $GFF3, '<', $gff3;

while ( my $input = <$GFF3> ) {
    chomp $input;
    if ( $input =~ /\A [#]/xms ) {
        print "$input\n";
    }
    elsif ( ( $input =~ /\A [^\t]* (?: \t [^\t]*){7} \t ID=([^\s;]+); /xms ) and ( exists $ok_genes{$1} ) ) {
        print "$input\n";
    }
    elsif ( ( $input =~ /\A [^\t]* (?: \t [^\t]*){7} \t .* gene_id=([^\s;]+); /xms ) and ( exists $ok_genes{$1} ) ) {
        print "$input\n";
    }
    elsif ( ( $input =~ /\A [^\t]* (?: \t [^\t]*){7} \t .* Parent=([^\s;]+); /xms ) and ( exists $ok_genes{$1} ) ) {
        print "$input\n";
    }
}
close $GFF3;

