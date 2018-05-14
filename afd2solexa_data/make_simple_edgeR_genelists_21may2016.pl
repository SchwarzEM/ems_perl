#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;
use List::MoreUtils qw(uniq);
use Scalar::Util qw(looks_like_number);

my %edge2deseq = (
    'ColWT' => 'Col_WT',
    'atml1' => 'atml1-3',
    'lgo' => 'lgo-2',
    'LGOoe' => 'ATML1__LGO',
    'LGOoe.atml1' => 'ATML1__LGO_atml1-3',
    'batch1' => 'batch_one',
    'batch2' => 'batch_two',
    'batch3' => 'batch_three',
);

# Sample (small) input file contents of CSVs_20may2016/ColWT_atml1_edgeR_exactTest_padj0.1_2016.05.20.tsv.txt
# Gene	logFC	logCPM	PValue	FDR
# AT1G53480	-5.91404088512938	1.85739067790375	1.31094090493115e-09	2.44162743543426e-05
# AT1G70890	2.25879380094297	2.458833024098	1.12896557778915e-06	0.0105134919431615

while (my $infile = <>) {
    chomp $infile;

    if (! -r $infile) {
        die "Cannot read file: $infile\n";
    }

    my $basename = basename ($infile);

    my $top    = q{};
    my $bottom = q{};
    if ( $basename =~ /\A (\S+) _ (\S+) _edgeR_exactTest_padj0\.1_2016\.05\.20\.tsv\.txt \z/xms ) {
        $top    = $1;
        $bottom = $2;
        if ( (! exists $edge2deseq{$top} ) or (! exists $edge2deseq{$bottom} ) ) {
            die "Cannot convert condition name(s) in: $infile\n";
        }
    }
    else {
        die "Cannot parse basename $basename\n";
    }    

    my $tag    = $edge2deseq{$top} . q{.vs.} . $edge2deseq{$bottom};

    my $up_outfile = "$tag.up.edgeR.gene_list.txt";
    $up_outfile = safename($up_outfile);

    my $down_outfile = "$tag.down.edgeR.gene_list.txt";
    $down_outfile = safename($down_outfile);

    my @up_tair_genes   = ();
    my @down_tair_genes = ();

    open my $INFILE, '<', $infile;

    while (my $input = <$INFILE>) {
        chomp $input;
        if ( $input =~ /\A (AT (?: 1|2|3|4|5|C|M) G\d+) \t (\S+) \t \S+ \t \S+ \t \S+ \z/xms ) {
            my $tair_gene       = $1;
            my $log2fold_change = $2;
            if (! looks_like_number($log2fold_change) ) {
                die "In input file $infile, cannot parse putative log2FoldChange \"$log2fold_change\" in: $input\n";
            }
            if ( $log2fold_change > 0 ) {
                push @up_tair_genes, $tair_gene;
            }
            if ( $log2fold_change < 0 ) {
                push @down_tair_genes, $tair_gene;
            }
            if ( $log2fold_change == 0 ) {
                die "Mars is invading, prepare to die!\n";
            }
        }
        elsif ( $input !~ /\A Gene \t logFC \t logCPM \t PValue \t FDR \z/xms ) {
            warn "In input file $infile, choose not to parse: $input\n";
        }
    }
    close $INFILE;

    @up_tair_genes = sort @up_tair_genes;
    @up_tair_genes = uniq @up_tair_genes;

    @down_tair_genes = sort @down_tair_genes;
    @down_tair_genes = uniq @down_tair_genes;

    open my $UP, '>', $up_outfile;
    foreach my $up_tair_gene (@up_tair_genes) {
        print $UP "$up_tair_gene\n";
    }
    close $UP;

    open my $DOWN, '>', $down_outfile;  
    foreach my $down_tair_gene (@down_tair_genes) {
        print $DOWN "$down_tair_gene\n";
    }
    close $DOWN;
}

sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}

