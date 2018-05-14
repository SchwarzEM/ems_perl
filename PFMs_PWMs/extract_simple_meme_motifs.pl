#!/usr/bin/env perl

# extract_simple_meme_motifs.pl -- Erich Schwarz <emsch@caltech.edu>, 11/10/2011.
# Purpose: given a bunch of ASCII text with MEME-style text motifs buried in it, print only the motifs.

use strict;
use warnings;

my $print_white_space = 1;
my $print_bg_vals     = 0;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A MEME [ ] version [ ] \d\S* .* \z /xms ) {
        print "$input\n";
        $print_white_space++;
    }
    elsif ( $input =~ /\A ALPHABET= [ ] ACGT \s* \z /xms ) { 
        print "$input\n";
        $print_white_space++;
    }
    elsif ( $input =~ /\A strands: [ ] \+ [ ] \- \s* \z/xms ) { 
        print "$input\n";
        $print_white_space++;
    }
    elsif ( $input =~ /\A Background [ ] letter [ ] frequencies [ ] \( from [ ] .+ \) : \s* \z/xms ) {
        print "$input\n";
        $print_white_space++;
        $print_bg_vals++;
    }
    elsif ( $input =~ /\A \s* A \s \d\.\d+ \s C \s \d\.\d+ \s  G \s \d\.\d+ \s T \s \d\.\d+ \s* \z/xms ) {
        print "$input\n" if $print_bg_vals;
        $print_white_space++;
        $print_bg_vals = 0;
    }
    elsif ( $input =~ /\A MOTIF [ ]+ \S .* \z/xms ) { 
        print "$input\n";
        $print_white_space++;
    }
    elsif ( $input =~ /\A BL \s+ MOTIF \s+ \S+ \s+ width=0 \s+ seqs=0 \s* \z/xms ) { 
        print "$input\n";
        $print_white_space++;
    }
    elsif ( $input =~ /\A letter-probability [ ] matrix: [ ] alength= [ ] 4 [ ] w= [ ] \d+ [ ] nsites= [ ] \d+ [ ] E= [ ] \S+ \s* \z/xms ) {
        print "$input\n";
        $print_white_space++;
    }
    elsif ( $input =~ /\A \s* \d\.\d+ \s+ \d\.\d+ \s+ \d\.\d+ \s+ \d\.\d+  \s* \z/xms ) {
        print "$input\n";
        $print_white_space++;
    }
    elsif ( $input =~ /\A \s* \z/xms ) {
        if ( $print_white_space ) {  
            print "$input\n";
            $print_white_space = 0;
        }
    }
}

