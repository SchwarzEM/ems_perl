#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (.+) \z/xms ) {
        my $gene_orig = $1;
        my $text      = $2;
        my $gene      = $gene_orig;

        if ( exists $data_ref->{$gene_orig} ) {
            my $i = $data_ref->{$gene_orig};
            $i++;
            $data_ref->{$gene_orig} = $i;
            $gene = $gene_orig . '_' . $i;
        }
        else {
            $data_ref->{$gene_orig} = 1;
        }
 
        print "$gene\t$text\n";
    }
    else {
    }
}

