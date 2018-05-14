#!/usr/bin/env perl

# patch_meme_motif_BL_line.pl -- Erich Schwarz <emsch@caltech.edu>, 11/10/2011.
# Purpose: patch a defect which seems to occur in some MEME motif texts.

use strict;
use warnings;

my $motif = q{};

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A MOTIF \s+ (\S+) .* /xms ) { 
        $motif = $1;
        print "$input\n";
        print "BL   MOTIF $motif width=0 seqs=0\n";
    }
    else { 
        print "$input\n";
    }
}

