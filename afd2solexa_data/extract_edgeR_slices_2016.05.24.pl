#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $label   = q{};
my $outfile = q{};
my @infiles = @ARGV;

chomp @infiles;

foreach my $infile (@infiles) {
    if ( $infile =~ /\A ( (\S+) _ (\S+) _edgeR_exactTest_padj0\.1_2016\.05\.20) \.tsv\.txt \z/xms ) { 
        $outfile    = $1;
        my $top     = $2;
        my $bottom  = $3;

        $label      = "$top.vs.$bottom";

        $outfile = $outfile . '.subtable.tsv.txt';
        $outfile = safename($outfile); 
    } 
    else {
        die "Cannot parse input file name: $infile\n";
    }
    open my $INFILE, '<', $infile;
    open my $OUTFILE, '>', $outfile;
    while (my $input = <$INFILE>) {
        chomp $input;

        $input =~ s/\A["]{2}/Gene/;
        $input =~ s/["]//g;
        $input =~ s/logFC/log2FoldChange[$label]/;
        $input =~ s/FDR/padj[$label]/;

        if ($input =~ /\A (\S+) \t (\S+) (?: \t \S+){2} \t (\S+) \z/xms ) {
            my $gene = $1;
            my $fold = $2;
            my $padj = $3;
            print $OUTFILE "$gene\t$fold\t$padj\n";
        }
        else {
            die "From input file $infile, cannot parse line: $input\n";
        }
    }
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

