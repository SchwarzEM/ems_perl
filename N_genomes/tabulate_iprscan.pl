#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $gene2cds = q{};
$gene2cds    = shift @ARGV if @ARGV;

my @infiles = ();
@infiles    = @ARGV if @ARGV;

if ( (! $gene2cds) or (! @infiles) ) {
    die "Format: tabulate_iprscan.pl [gene-to-CDS.tsv.txt] [1+ InterProScan *.tsv outputs]\n";
}

my @genes  = ();
my $header = "Gene\tInterPro\tGO_term";

my $data_ref;

open my $GENE2CDS, '<', $gene2cds;
while (my $input = <$GENE2CDS>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $gene       = $1;
        my $transcript = $2;
        push @genes, $gene;
        $data_ref->{'tx'}->{$transcript}->{'gene'} = $gene;
    }
}
close $GENE2CDS;

foreach my $infile (@infiles) {
    my $gene       = q{};
    my $transcript = q{};

    open my $INFILE, '<', $infile;
    while (my $input = <$INFILE>) { 
        chomp $input;
        if ( $input =~ /\A (\S+) \t (?: [^\t]* \t){10} (IPR\d+) \t ([^\t]+) \t (\S.*\S) \b /xms ) { 
            $transcript  = $1;
            my $ipr_acc  = $2;
            my $ipr_desc = $3;
            my $go_text  = $4;
            $gene = $data_ref->{'tx'}->{$transcript}->{'gene'};

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
    close $INFILE;
}

@genes = uniq(@genes);

foreach my $gene (@genes) { 
    my @ipr_motifs     = ();
    my $ipr_motif_text = q{};
    if ( exists $data_ref->{'gene'}->{$gene}->{'interpro'} ) { 
        @ipr_motifs = sort keys %{ $data_ref->{'gene'}->{$gene}->{'interpro'} };
        $ipr_motif_text = join '; ', @ipr_motifs;
    }

    my @go_terms     = ();
    my $go_term_text = q{};
    if ( exists $data_ref->{'gene'}->{$gene}->{'go_term'} ) {
        @go_terms = sort keys %{ $data_ref->{'gene'}->{$gene}->{'go_term'} };
        $go_term_text = join '; ', @go_terms;
    }

    print "$header\n" if $header;
    $header = q{};
    print "$gene\t$ipr_motif_text\t$go_term_text\n";
}

