#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $infile = q{};
$infile    = $ARGV[0] if $ARGV[0];

if (! $infile ) {
    die "Format: split_annovar_annots_11mar2026.pl [input ANNOVAR mutation annots] > [split between shared genes with added 'shared' caveat]\n";
}

open my $INFILE, '<', $infile;

while ( my $input = <$INFILE> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (.+) \z/xms ) {
        my $orig_gene_id = $1;
        my $orig_annot   = $2;
        if ( $orig_gene_id !~ /[,]/xms ) { 
            print "$orig_gene_id\t$orig_annot\n";
        }
        else {
            my @gene_ids = split /[,]/, $orig_gene_id;
            foreach my $gene_id (@gene_ids) {
                print "$gene_id\t$orig_annot, shared between $orig_gene_id\n";
            }
        }
    }
}
close $INFILE;


