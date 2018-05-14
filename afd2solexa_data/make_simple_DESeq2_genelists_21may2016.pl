#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;
use List::MoreUtils qw(uniq);
use Scalar::Util qw(looks_like_number);

# Sample (small) input file contents:
# Gene  log2FoldChange[Col_WT.vs.atml1-3]       padj[Col_WT.vs.atml1-3]
# AT1G53480     -2.04234405137481       9.51320542899225e-09
# AT1G70890     1.36176060639264        0.0227251484079386

while (my $infile = <>) {
    chomp $infile;

    if (! -r $infile) {
        die "Cannot read file: $infile\n";
    }

    my $basename = basename ($infile);

    my $tag = q{};
    if ( $basename =~ /\A (\S+) _DESeq2_q0\.1_3nt_2016\.02\.13.subtable\.tsv\.txt \z/xms ) { 
        $tag = $1;
    }
    else {
        die "Cannot parse basename $basename\n";
    }    

    my $up_outfile = "$tag.up.gene_list.txt";
    $up_outfile = safename($up_outfile);

    my $down_outfile = "$tag.down.gene_list.txt";
    $down_outfile = safename($down_outfile);

    my @up_tair_genes   = ();
    my @down_tair_genes = ();

    open my $INFILE, '<', $infile;

    while (my $input = <$INFILE>) {
        chomp $input;
        if ( $input =~ /\A (AT (?: 1|2|3|4|5|C|M) G\d+) \t (\S+) \t \S+ \z/xms ) {
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
        elsif ( $input !~ /\A Gene \t log2FoldChange\[ \S+ \] \t padj\[ \S+ \] \z/xms ) {
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

