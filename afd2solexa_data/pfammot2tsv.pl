#!/usr/bin/env perl

# pfammot2tsv.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/23/2010.
# Purpose: given a Tablemaker output of PFAM (etc.) motifs appended to their genes, make a .tsv with one line per gene and PFAM only.

use strict;
use warnings;

my $gene             = q{};
my $motif_id         = q{};
my $motif_desc       = q{};
my $genes2motifs_ref;

# Sample input lines:
# 
# "WBGene00000006"        "INTERPRO:IPR004841"    "Amino acid permease domain"
# "WBGene00000006"        "PFAM:PF00324"  "Amino acid permease"
# "WBGene00000006"        "PhosphoPep:C55C2.5a"   "Phosphorylation site"

while (my $input = <>) { 
    chomp $input;
    # To get both types of motif, use '((?:INTERPRO|PFAM):\w+)':
    if ( $input =~ /\A \" (WBGene\d+) \" \s+ \"PFAM:(\w+)\" \s+ \" ([^\"]+) \" \s* \z /xms ) { 
        $gene       = $1;
        $motif_id   = $2;
        $motif_desc = $3;
        # Get rid of obnoxious escape-characterization of forward slashes (meant to deal with ACeDB!):
        $motif_desc =~ s/\\//g;
        my $text = $motif_id . ': ' . $motif_desc;
        $genes2motifs_ref->{$gene}->{$text} = 1;
    }
}

foreach my $g1 ( sort keys %{ $genes2motifs_ref } ) { 
    my @annots = sort keys %{ $genes2motifs_ref->{$g1} };
    my $annot_line = join '; ', @annots;
    print "$g1\t\"$annot_line\"\n";
}

