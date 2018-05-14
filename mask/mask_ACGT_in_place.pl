#!/usr/bin/env perl

# mask_ACGT_in_place.pl -- Erich Schwarz, 10/8/10, <emsch@its.caltech.edu>.
# Purpose: use 'N' to mask uppercase nucleotides in a FASTA file (keeping the filename, while masking its contents).

use strict;
use warnings;

local $^I = '.bak';

while (my $input = <>) {
    chomp $input;
    unless ($input =~ /^>/) {
        $input =~ s/[A-Z]/N/g;
    }
    print "$input\n";
}

