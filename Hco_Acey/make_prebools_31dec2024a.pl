#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $i = 0;

while ( my $input = <> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \.clusts\.\S+\.tsv.txt/xms ) {
        my $tag = $1;
        my $outfile = 'Sherm.Pfam-36.0.vs.' . $tag . '.tsv.txt';
        $i++;
        print 'add_tab_annots.pl -f ',
              '-i ../annots/Sherm_2024.09.30.genes.txt ',
              '-a ../annots_prev/sherm_braker_combined_2022.12.17.hmmscan_cut_ga.Pfam-36.0.tsv.txt ',
              "$input ",
              '1>',
              "$outfile ",
              '2>',
              "test$i.err ;",
              "\n",
              ;
    }
    else {
        die "Cannot parse input: $input\n";
    }
}
