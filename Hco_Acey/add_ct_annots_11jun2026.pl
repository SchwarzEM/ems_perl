#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %ct2txts     = ();
my $ct_annots   = q{};
my $gene_annots = q{};

$ct_annots   = $ARGV[0] if $ARGV[0];
$gene_annots = $ARGV[1] if $ARGV[1];

if ( (! $ct_annots ) or (! $gene_annots ) ) {
    die "Format: add_ct_annots_11jun2026.pl [CT annots] [gene/orth annots] > [gene/CT annots]\n";
}

open my $CT, '<', $ct_annots;
while ( my $input = <$CT> ) {
    chomp $input;
    my $rem_id = q{};
    my $ct_txt = q{};
    if ( $input =~ /\A (\S+) \t ([^\t]+) \z/xms ) {
        $rem_id = $1;
        $ct_txt = $2;
        $ct2txts{$rem_id} = "$ct_txt [$rem_id]";
    }
    else {
        die "From CT annots $ct_annots, cannot parse: $input\n";
    }
}
close $CT;

my @rem_ids = sort keys %ct2txts;

open my $GENE, '<', $gene_annots;
while ( my $input = <$GENE> ) {
    chomp $input; 
    my $gene    = q{};
    my $orth    = q{};
    my %ct_desc = ();
    if ( $input =~ /\A (\S+) \t ([^\t]+) \z/xms ) {
        $gene = $1;
        $orth = $2;
        foreach my $rem_id (@rem_ids) {
            if ( $orth =~ /$rem_id/ ) {
                my $ct_txt = $ct2txts{$rem_id};
                $ct_desc{$ct_txt} = 1;
            }
        }
        my @descs = sort keys %ct_desc;
        my $desc_txt = join '; ', @descs;
        if ( $gene eq 'Gene' ) {
            print "Gene\tThomas_2010_annot\n";
        }
        else {
            print "$gene\t$desc_txt\n";
        }
    }
    else {
        die "From gene annots $gene_annots, cannot parse: $input\n";
    }
}
close $GENE;

