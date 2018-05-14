#!/usr/bin/env perl

# uniq_GFF_seqlines.pl -- Erich Schwarz <emsch@its.caltech.edu>, 7/11/2008.
# Purpose: have one sequence line per segment in GFF; act as filter.

use strict;
use warnings;

my $sequence     = q{};
my %stored_lines = ();
my %do_not_store = ();

while (my $input = <>) { 
    chomp $input;

    # Store ## lines w/o immediate printing.
    # E.g.: ##sequence-region Contig1000 1 15788
    if ( $input =~ / \A 
                        [#]{2}
                        sequence-region 
                        \s+ 
                        (.+?) 
                        \s+ 
                        1 
                        \s 
                        \d+ /xms ) { 
        print "$input\n";
    }

    elsif ( $input =~ / \A 
                        ([^\t]+) 
                        \t 
                        \. 
                        \t 
                        Sequence 
                        \t 
                        1 
                        \t 
                        \d+ /xms ) { 
        print "$input\n";
    }
}

