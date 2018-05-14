#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use List::MoreUtils qw(uniq);

my $gene_ids = $ARGV[0];
my $data     = $ARGV[1];

my $header   = "Gene\tKOG/TWOG/LSE\tKOG_code\n";
my @output   = ();

# Note that we *have* to allow redundant mapping of CDS IDs to gene IDs, because we had the peptides-are-non-unique problem.
my $cds2gene_ref = ();

open my $GENE_IDS, '<', $gene_ids;
while (my $input = <$GENE_IDS>) {
   chomp $input;

# Sample input:
# WBGene00000001,aap-1,Y110A7A.10
# WBGene00000002,aat-1,F27C8.1
# WBGene00000087,aex-4,
# WBGene00000230,,
# WBGene00000343,,G44K07.8

    if ( $input =~ /\A (WBGene\d+) [,] (.*) [,] ([^,\s]+) /xms ) { 
        my $wb_id  = $1;
        my $cgc_id = $2;
        my $cds_id = $3;
        my $gene_id = "$wb_id|$cds_id";
        if ($cgc_id) {
            $gene_id .= "|$cgc_id";
        }
        $cds2gene_ref->{$cds_id}->{$gene_id} = 1;
    }
}
close $GENE_IDS;

open my $DATA, '<', $data;
while (my $input = <$DATA>) {
    chomp $input;

# Sample input:
# 2L52.1	KOG1721 [Zn-finger]	Code_R
# 6R55.1	KOG3989 [Beta-2-glycoprotein I]	Code_W
# AC3.4	KOG1315 [Predicted DHHC-type Zn-finger protein]; LSE0126 [Uncharacterized protein]	Code_R; Code_S

    if ( $input =~ /\A (\S+) \t (.+) \z/xms ) { 
        my $cds_id = $1;
        my $annot  = $2;
        if (! exists $cds2gene_ref->{$cds_id} ) {
            die "Cannot map CDS $cds_id to full gene name in WS130\n";
        }
        my @gene_ids = sort keys %{ $cds2gene_ref->{$cds_id} };
        foreach my $gene_id (@gene_ids) {
            my $output_line = "$gene_id\t$annot\n";
            push @output, $output_line;
        }
    }
}
close $DATA;

@output = sort @output;
@output = uniq @output;
print $header;
print @output;


