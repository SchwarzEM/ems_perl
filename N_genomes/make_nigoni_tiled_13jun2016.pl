#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $tmp_dir = 'tmp';
$tmp_dir = safename($tmp_dir);

my $header1 = '#!/bin/bash';
my $header2 = "    mkdir $tmp_dir ;";

my $infile   = $ARGV[0];
my $out_desc = $ARGV[1];

if ( $out_desc =~ /\W/xms ) { 
     die "Cannot use output description \"$out_desc\"\n";
}

open my $INFILE, '<', $infile;
while (my $input = <$INFILE>) {
    chomp $input;
    if ( $input =~ / \A (nigoni_2015.12.01_\S+) \t ([+]|[-]) \z/xms ) { 
        my $contig = $1;
        my $sense  = $2;

        print "$header1\n\n" if $header1;
        $header1 = q{};

        print "$header2\n\n" if $header2;
        $header2 = q{};

        print "    extract_fasta_subset.pl -f nigoni_2015.12.01.rev1_gen.dna.fa -l $contig > tmp/$contig.fa ;\n";

        if ( $sense eq '+' ) {
            print "    cat tmp/$contig.fa >> tmp/nigoni_2015.12.01.rev1_tiled.$out_desc.contigs_gen.dna.fa ;\n";
        }
        elsif ( $sense eq '-' ) {
            print "    make_revcomp_seqs.pl -s -f tmp/$contig.fa | tag_FASTA_names.pl -i - -s \"_revcomp\" > tmp/$contig.rc ;\n";
            print "    cat tmp/$contig.rc >> tmp/nigoni_2015.12.01.rev1_tiled.$out_desc.contigs_gen.dna.fa ;\n";
        }

        print "\n";
    }
    else {
        die "Cannot parse: $input\n";
    }
}
close $INFILE;

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

