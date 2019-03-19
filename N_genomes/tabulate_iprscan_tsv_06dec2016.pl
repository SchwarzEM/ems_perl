#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my $gene = q{};

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \.t\d+ \t (?: [^\t]* \t){10} (IPR\d+) \t ([^\t]+) \t (\S.*\S) \b /xms ) { 
        my $gene     = $1;
        my $ipr_acc  = $2;
        my $ipr_desc = $3;
        my $go_text  = $4;

        my $full_ipr_desc = "$ipr_desc [$ipr_acc]";
        $data_ref->{'gene'}->{$gene}->{'interpro'}->{$full_ipr_desc} = 1;

        if ( $go_text =~ /\A GO[:]\d+ /xms ) { 
            my @go_terms = split /\|/, $go_text;
            foreach my $go_term (@go_terms) {
                $data_ref->{'gene'}->{$gene}->{'go_term'}->{$go_term} = 1;
            }
        }
    }
}

my $header = "Gene\tInterPro\tGO_term";

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene1 (@genes) { 
    my @ipr_motifs     = ();
    my $ipr_motif_text = q{};
    if ( exists $data_ref->{'gene'}->{$gene1}->{'interpro'} ) { 
        @ipr_motifs = sort keys %{ $data_ref->{'gene'}->{$gene1}->{'interpro'} };
        $ipr_motif_text = join '; ', @ipr_motifs;
    }

    my @go_terms     = ();
    my $go_term_text = q{};
    if ( exists $data_ref->{'gene'}->{$gene1}->{'go_term'} ) {
        @go_terms = sort keys %{ $data_ref->{'gene'}->{$gene1}->{'go_term'} };
        $go_term_text = join '; ', @go_terms;
    }

    print "$header\n" if $header;
    $header = q{};
    print "$gene1\t$ipr_motif_text\t$go_term_text\n";
}

