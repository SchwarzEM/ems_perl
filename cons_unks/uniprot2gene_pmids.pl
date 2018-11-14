#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my $uni2pmid = q{};
my $uni2gene = q{};

$uni2gene = $ARGV[0] if	$ARGV[0];
$uni2pmid = $ARGV[1] if	$ARGV[1];

if ( (! $uni2pmid) or (! $uni2gene) ) {
    die "Format: uniprot2gene_pmids.pl [uniprot2gene] [uniprot2pmid] > [genes, PMIDs (and weights), summed weighted annotation]\n";
}

open my $UNI2GENE, '<', $uni2gene;
while (my $input = <$UNI2GENE>) {
    chomp $input;
    # For some weird reason, I can run this with (.+), but not with ([^\t]+), which is what I really wanted.  See below for work-around.
    if ( $input =~ /\A (\S+) \t (.+) \z/xms ) {
        my $uniprot   = $1;
        my $gene_text = $2;

        # This is a work-around to deal with my being unable to screen directly for [^\t]+, above.
        if ( $gene_text =~ /\t/xms ) {
            warn "Unwanted tab? inside \"$gene_text\", in: $input\n";
        }

        my @genes     = split /; /, $gene_text;
        foreach my $gene (@genes) {
            $data_ref->{'uniprot'}->{$uniprot}->{'gene'}->{$gene} = 1;
        }
    }
    else {
	die "In uni2gene file $uni2gene, cannot parse: $input\n";
    }
}
close $UNI2GENE;

open my $UNI2PMID, '<', $uni2pmid;
while (my $input = <$UNI2PMID>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t ([^t]+) \z/xms ) {
        my $uniprot   = $1;
        my $pmid_text = $2;
        my @pmids = split /; /, $pmid_text;

        if ( exists $data_ref->{'uniprot'}->{$uniprot}->{'gene'} ) {
            my @genes = sort keys %{ $data_ref->{'uniprot'}->{$uniprot}->{'gene'} };
            foreach my $gene (@genes) {
                foreach my $pmid (@pmids) {
                    # gene to pmid, obviously
                    $data_ref->{'gene'}->{$gene}->{'pmid'}->{$pmid} = 1;
                    # but also: get a count of each PMID's citation count within a given proteome
                    $data_ref->{'pmid'}->{$pmid}->{'citations'}++;
                }
            }
        }

        else {
            warn "Could not map uniprot $uniprot to a gene\n";
        }
    }
    else {
        die "In uni2pmid file $uni2pmid, cannot parse: $input\n";
    }
}
close $UNI2PMID;

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene (@genes) {
    my @pmids           = sort keys %{ $data_ref->{'gene'}->{$gene}->{'pmid'} };
    my $annot_sum       = 0;
    my @pmids_w_weights = ();
    my $pmid_text       = q{};

    foreach my $pmid (@pmids) {
        my $citations      = $data_ref->{'pmid'}->{$pmid}->{'citations'};
        my $weighted_annot = (1/$citations);
        my $pmid_w_weight  = "$pmid [$weighted_annot]";

        $annot_sum = $annot_sum + $weighted_annot;
        push @pmids_w_weights, $pmid_w_weight;
    }
    @pmids_w_weights = sort @pmids_w_weights;
    $pmid_text       = join '; ', @pmids_w_weights;     
    $data_ref->{'final_gene'}->{$gene}->{'annot_sum'} = $annot_sum;
    $data_ref->{'final_gene'}->{$gene}->{'pmid_text'} = $pmid_text;
}

# Do Schwartzian transformation to get the final table sorted by descending $annot_sum values.

my @final_genes = sort {     $data_ref->{'final_gene'}->{$b}->{'annot_sum'} 
                         <=> $data_ref->{'final_gene'}->{$a}->{'annot_sum'} } 
                       keys %{ $data_ref->{'final_gene'} };

foreach my $gene (@final_genes) {
    my $pmid_text = $data_ref->{'final_gene'}->{$gene}->{'pmid_text'};
    my $annot_sum = $data_ref->{'final_gene'}->{$gene}->{'annot_sum'};
    print "$gene\t$pmid_text\t$annot_sum\n";
}

