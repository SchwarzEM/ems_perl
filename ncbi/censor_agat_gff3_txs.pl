#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $gff3     = q{};
my $ban_list = q{};
my %banned   = ();

$gff3     = $ARGV[0] if $ARGV[0];
$ban_list = $ARGV[1] if $ARGV[1];

if ( (! $gff3 ) or (! $ban_list ) ) {
    die "Format: censor_agat_gff3_txs.pl [GFF3] [banned_tx_list] > [filtered GFF3]\n";
}

open my $BAN, '<', $ban_list;
while ( my $tx = <$BAN> ) {
    chomp $tx;
    $banned{$tx} = 1;
}
close $BAN;

open my $GFF3, '<', $gff3;
while ( my $annot = <$GFF3> ) {
    chomp $annot;

    # [Sample input lines:]
    # Sherm_X AUGUSTUS        transcript      8236    8835    0.29    -       .       ID=Sherm_X.g16153.t1;Parent=Sherm_X.g16153
    # [or:]
    # Sherm_X AUGUSTUS        exon    8236    8835    .       -       .       ID=exon-217435;Parent=Sherm_X.g16153.t1;gene_id=Sherm_X.g16153;transcript_id=Sherm_X.g16153.t1
    # Sherm_X AUGUSTUS        CDS     8236    8835    0.29    -       0       ID=cds-217433;Parent=Sherm_X.g16153.t1;gene_id=Sherm_X.g16153;transcript_id=Sherm_X.g16153.t1
    # Sherm_X AUGUSTUS        start_codon     8833    8835    .       -       0       ID=start_codon-34914;Parent=Sherm_X.g16153.t1;gene_id=Sherm_X.g16153;transcript_id=Sherm_X.g16153.t1
    # Sherm_X AUGUSTUS        stop_codon      8236    8238    .       -       0       ID=stop_codon-34917;Parent=Sherm_X.g16153.t1;gene_id=Sherm_X.g16153;transcript_id=Sherm_X.g16153.t1

    if ( $annot =~ /\A \S+ \t \S+ \t transcript \t .+ \t ID= ([^;\s]+)/xms ) {
        my $tx = $1;
        if (! exists $banned{$tx} ) {
            print "$annot\n";
        }
    }
    elsif ( $annot =~ /\A .+ \t ID=[^\t]+ transcript_id= ([^;\s]+)/xms ) {
        my $tx = $1;
        if (! exists $banned{$tx} ) {
            print "$annot\n";
        }
    }
    else {
        print "$annot\n";
    }
}
close $GFF3;

