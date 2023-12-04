#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $infile = q{};
my $patch  = q{};

$infile = $ARGV[0] if $ARGV[0];
$patch  = $ARGV[1] if $ARGV[1];

if ( (! $infile ) or (! $patch ) ) {
    die "Format: patch_second_fisher_motif_03dec2023.pl [infile] [2cd motif to patch] > [outfile]\n";
}

open my $INFILE, '<', $infile;
while ( my $input = <$INFILE> ) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \t (.*) \z/xms ) {
        my $motif = $1;
        my $text  = $2;
        if ( $motif eq 'Motif' ) {
            $text =~ s/Motif_genes/Motif_1_genes/g;
            $text =~ s/Class_genes/Motif_2_genes/g;
            $text =~ s/Motif.Class_overlap/Motif_1.Motif_2_overlap/g;
            print "Motif_1\tMotif_2\t$text\n";
        }
        else {
            print "$motif\t$patch\t$text\n";
        }
    }
    else {
        die "From input file $infile, cannot format: $input\n";
    }
}
close $INFILE;
