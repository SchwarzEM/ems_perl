#!/usr/bin/env perl

# get_unmap.genes_w_flags.pl -- Erich Schwarz <ems394@cornell.edu>, 6/18/2018.
# Purpose: given *.gtf.unmap from liftover, extract only 'gene' lines with immediately preceding comment lines.

use strict;
use warnings;
use autodie;

my $comment = q{};
my $text    = q{};

while (my $input = <>) {
    chomp $input;
    $comment = $text;
    $text    = $input;
    if ( $input =~ /\A [^\t]* \t [^\t]* \t gene \t /xms ) { 
        print "$comment\n";
        print "$text\n";
    }
}


