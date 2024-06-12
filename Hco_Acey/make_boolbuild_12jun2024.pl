#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $i = 0;
my $j = q{};

while (my $input = <>) {
    chomp $input;
    # Sample input: Nippo_ES.L4_annots_2024.05.09.01.tsv.txt
    if ( $input =~ /\A (\S+) _annots .* \.tsv\.txt \z/xms ) {
        my $stem = $1;

        $i++;
        $j = sprintf "%02i", $i;

        print "\n";
        print 'add_tab_annots.pl -f -i ../annots/nippo_braker_combined_2023.07.15.genes.txt';
        print ' -a ../annots/nippo_braker_combined_2023.07.15.pep.hmmscan_cut_ga.Pfam-36.0.tbl.tsv.txt';
        print " $input";
        print " 1>Nippo.Pfam-36.0.vs.$stem.tsv.txt 2>test$j.err ;";
        print "\n\n";

        $i++;
        $j = sprintf "%02i", $i;

        print 'add_tab_annots.pl -f -i ../annots/nippo_braker_combined_2023.07.15.genes.txt';
        print ' -a ../annots/Nippo.2023.07.15.interproscan-5.67-99.0.tsv.txt';
        print " $input";
        print " 1>Nippo.interproscan-5.67-99.0.vs.$stem.tsv.txt 2>test$j.err ;";
        print "\n\n";

        $i++;
        $j = sprintf "%02i", $i;

        print 'add_tab_annots.pl -f -i ../annots/nippo_braker_combined_2023.07.15.genes.txt';
        print ' -a ../annots/nippo_2023.07.15.phobius.short.tsv.txt';
        print " $input";
        print " 1>Nippo.phobius.vs.$stem.tsv.txt 2>test$j.err ;";
        print "\n";
    }
    else {
        die "Cannot parse input: $input\n";
    }
}

