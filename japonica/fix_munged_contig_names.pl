#!/usr/bin/perl

# fix_munged_contig_names.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/17/2008.
# Purpose: fix mis-named contigs with '.fa' suffixes that are jamming GFF loading.

use strict;
use warnings;

while (my $input = <>) { 
    if ($input =~ / \A ((Contig\d+)\.fa) /xms ) { 
        my $bad_name = $1;
        my $good_name = $2;
        $input =~ s/$bad_name/$good_name/g;
    }
    print $input;
}

