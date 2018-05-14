#!/usr/bin/env perl

# fasta_Xs.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/8/2009.
# Purpose: remove terminal '*' from seqs.; replace internal '*'s or '-'s with 'X's.

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if (    ( $input =~ /\A > /xms       ) 
         or ( $input =~ / \A \s* \z /xms ) ) { 
        print "$input\n";
    }
    else { 
        $input =~ s/\*\s*\z//;
        $input =~ s/\*/X/g;
        $input =~ s/\-/X/g;
        print "$input\n";
    }
}

