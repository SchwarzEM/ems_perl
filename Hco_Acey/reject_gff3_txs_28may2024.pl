#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $txlist = q{};
my $gff3     = q{};

my %bad_txs = ();

$txlist = $ARGV[0] if $ARGV[0];
$gff3     = $ARGV[1] if $ARGV[1];

if ( (! $txlist ) or (! $gff3 ) ) {
    die "Format: select_gff3_genes_28may2024.pl [transcript/mRNA/CDS list] [GFF3] > [GFF3 *without* listed txs/mRNAs/CDSes]\n";
}

open my $TXLIST, '<', $txlist;
while ( my $tx = <$TXLIST> ) {
    chomp $tx;
    $bad_txs{$tx} = 1;
}
close $TXLIST;

open my $GFF3, '<', $gff3;

while ( my $input = <$GFF3> ) {
    chomp $input;
    my $print = 1;
    if ( ( $input =~ /\A [^\t]* (?: \t [^\t]*){7} \t ID=([^\s;]+); /xms ) and ( exists $bad_txs{$1} ) ) {
        $print = 0;
    }
    elsif ( ( $input =~ /\A [^\t]* (?: \t [^\t]*){7} \t .* transcript_id=([^\s;]+); /xms ) and ( exists $bad_txs{$1} ) ) {
        $print = 0;
    }
    elsif ( ( $input =~ /\A [^\t]* (?: \t [^\t]*){7} \t .* Parent=([^\s;]+); /xms ) and ( exists $bad_txs{$1} ) ) {
        $print = 0;
    }
    if ( $print ) {
        print "$input\n";
    }
}
close $GFF3;

