#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Getopt::Long;

my @infiles = ();
my $gene2tx     = q{};
my $gene        = q{};

my $data_ref;

my $help;

GetOptions ( 
    'iprscan_tsv=s{,}' => \@infiles,
    'gene2tx=s'        => \$gene2tx,
    'help'             => \$help, 
);

if ( $help or (! -r $gene2tx ) or (! @infiles ) ) {
    die "tabulate_iprscan_tsv.pl\n",
        "    --iprscan_tsv|-i  [one or more TSV outputs from InterProScan]\n",
        "    --gene2tx|-g      [gene-to-transcript TSV table]\n",
        "    --help|-h         [print this message]\n",
        ;
}

open my $GENE2TX, '<', $gene2tx;
while (my $input = <$GENE2TX>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $gene = $1;
        my $tx   = $2;
        if ( exists $data_ref->{'tx'}->{$tx}->{'gene'} ) {
            die "Redundant mapping of transcript $tx to $data_ref->{'tx'}->{$tx}->{'gene'} and $gene\n";
        }
        $data_ref->{'tx'}->{$tx}->{'gene'} = $gene;
    }
    else { 
        die "From gene-to-transcript table $gene2tx, cannot parse: $input\n";
    }
}
close $GENE2TX;

foreach my $infile (@infiles) {
    open my $INFILE, '<', $infile;
    while (my $input = <$INFILE>) { 
        chomp $input;
        # Sample InterPro ID: IPR035500
        if ( $input =~ /\A (\S+) \t (?: [^\t]* \t){10} (IPR\d{6}) \t ([^\t]+) (?:\t|\z) /xms ) { 
            my $tx       = $1;
            my $ipr_acc  = $2;
            my $ipr_desc = $3;

            if (! exists $data_ref->{'tx'}->{$tx}->{'gene'} ) {
                die "Cannot map transcript $tx to a gene\n";
            }
            my $gene = $data_ref->{'tx'}->{$tx}->{'gene'};

            my $full_ipr_desc = "$ipr_desc [$ipr_acc]";
            $data_ref->{'gene'}->{$gene}->{'interpro'}->{$full_ipr_desc} = 1;
        }
        elsif ( $input =~ /IPR\d{6}/xms ) {
            die "Cannot parse input: $input\n";
        }
    }
    close $INFILE;
}

my $header = "Gene\tInterPro";

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene1 (@genes) { 
    my @ipr_motifs     = ();
    my $ipr_motif_text = q{};
    if ( exists $data_ref->{'gene'}->{$gene1}->{'interpro'} ) { 
        @ipr_motifs = sort keys %{ $data_ref->{'gene'}->{$gene1}->{'interpro'} };
        $ipr_motif_text = join '; ', @ipr_motifs;
    }
    print "$header\n" if $header;
    $header = q{};
    print "$gene1\t$ipr_motif_text\n";
}

