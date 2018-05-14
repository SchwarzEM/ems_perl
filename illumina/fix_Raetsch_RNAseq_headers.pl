#!/usr/bin/env perl

# fix_Raetsch_RNAseq_headers.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/6/2010.
# Purpose: fix amazingly badly named Raetsch RNAseq headers.

use strict;
use warnings;

my $header     = q{};
my $abs_suffix = 0;
my $abs_lines  = 0;
my $rel_suffix = 0;
my $rel_lines  = 0;
my $suffix     = 0;

while (my $input = <>) { 
    chomp $input;
    $abs_lines++;
    $rel_suffix = ($abs_suffix % 2);
    $rel_lines  = ($abs_lines % 4);

    if ( ( $input =~ /\A (\@ \d+) \s* \z/xms ) and ( $rel_lines == 1) ) {
        $header = $1;
        $abs_suffix++;
        $suffix = ($rel_suffix + 1);
        print $header, '_', "$suffix\n", ;
    }
    else { 
        print "$input\n";
    }
}

