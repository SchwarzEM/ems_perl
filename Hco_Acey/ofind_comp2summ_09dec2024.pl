#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $venom_list = q{};
my $ofinds     = q{};

$venom_list = $ARGV[0] if $ARGV[0];
$ofinds     = $ARGV[1] if $ARGV[1];

my %venom = ();

if ( (! $venom_list ) or (! $ofinds ) ) {
    die "Format: ofind_comp2summ_09dec2024.pl [venom gene list] [OrthoFinder complete descs.] > [summary with venoms noted]\n";
}

open my $VENOM_LIST, '<', $venom_list;
while ( my $input = <$VENOM_LIST> ) {
    chomp $input;
    if ( $input =~ /\A \S+ \z/xms ) {
        $venom{$input} = 1;
    }
    else {
        die "From venom list $venom_list, cannot parse: $input\n";
    }
}
close $VENOM_LIST;

open my $OFINDS, '<', $ofinds;

while ( my $input = <$OFINDS> ) {
    chomp $input;
    if ( $input =~ /\A ( O\S*\d+\(\d+[ ]genes,\d+[ ]taxa\)[:] \s+) (\S.+) \z /xms ) {
        my $ogroup_id_text = $1;
        my $gene_text      = $2;
        if ( $gene_text !~ /\S+\s*/xms ) {
            die "From OrthoFinder complete descriptions $ofinds, cannot parse gene text: $gene_text\n";
        }
        my %taxon_count = ();
        while ( $gene_text =~ /(\S+\s*)/xmsg ) {
            my $gene = $1;
            $gene =~ s/\s+\z//;
            if ( $gene =~ /\A (\S+) \( (\S+) \) \z/xms ) {
                my $gene_name = $1;
                my $taxon     = $2;
                if ( exists $venom{$gene_name} ) {
                    $taxon = "$taxon.venom";
                }
                $taxon_count{$taxon}++;
            }
            else {
                die "From OrthoFinder complete descriptions $ofinds, cannot parse gene: $gene\n";
            }
        }
        my @taxa = sort keys %taxon_count;
        my @descs = ();
        foreach my $taxon (@taxa) {
            my $desc = "$taxon ($taxon_count{$taxon} g.)";
            push @descs, $desc;
        }
        my $desc_text = join '; ', @descs;
        print "$ogroup_id_text$desc_text\n";
    }
    else {
        die "From OrthoFinder complete descriptions $ofinds, cannot parse: $input\n";
    }
}
close $OFINDS;

