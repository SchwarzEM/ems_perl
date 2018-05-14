#!/usr/bin/perl

# uncomment_FASTA_headers.pl -- Erich Schwarz, <emsch@its.caltech.edu>, 4/30/08
# Purpose: strip out all text other than '>Sequence_name' from FASTA headers.

use strict;
use warnings;
use Carp;

while (my $input = <>) { 
    chomp $input;
    if ($input !~ / \A > \S /xms) {
        print "$input\n";
    }
    elsif ($input =~ / \A > (\S+) .*? \z /xms) {
        my $input1 = $1;
        print '>';
        print $input1;
        print "\n";
    }
    else { 
        croak "Can't parse input line: $input\n";
    }
}

