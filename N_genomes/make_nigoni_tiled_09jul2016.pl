#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header1 = '#!/bin/bash';
my $header2 = '    mkdir tmp ;';

while (my $input = <>) {
    chomp $input;
    if ( $input =~ / \t ([+]|[-]) \t (nigoni_2016.07.08_\d+) \z/xms ) { 
        my $sense  = $1;
        my $contig = $2;

        print "$header1\n\n" if $header1;
        $header1 = q{};

        print "$header2\n\n" if $header2;
        $header2 = q{};

        print "    extract_fasta_subset.pl -f nigoni_2016.07.08_gen.dna.fa -l $contig > tmp/$contig.fa ;\n";

        if ( $sense eq '+' ) {
            print "    cat tmp/$contig.fa >> tmp/nigoni_2016.07.08.rev1_tiled.scaffolds_gen.dna.fa ;\n";
        }
        elsif ( $sense eq '-' ) {
            print "    make_revcomp_seqs.pl -s -f tmp/$contig.fa | tag_FASTA_names.pl -i - -s \"_revcomp\" > tmp/$contig.rc ;\n";
            print "    cat tmp/$contig.rc >> tmp/nigoni_2016.07.08.rev1_tiled.scaffolds_gen.dna.fa ;\n";
        }

        print "\n";
    }
}

