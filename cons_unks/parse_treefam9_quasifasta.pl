#!/usr/bin/env perl

# parse_treefam9_quasifasta.pl -- Erich Schwarz <ems394@cornell.edu>, 8/12/2013.
# Purpose: convert TreeFam 9.0's quasi-FASTA into something both somewhat correctly formatted and easy to extract subsets from.

use strict;
use warnings;

my $tf_id = q{};

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (TF\d+) /xms ) { 
        $tf_id = $1;
    }
    elsif ( $input =~ /\A > (\S.+) \z/xms ) {
        my $text = $1;
        print '>', $tf_id, '_', "$text\n";
    }
    elsif ( $input =~ /\A > /xms ) {
        die "Can't parse input line: $input\n";
    }
    elsif ( $input !~ /\A [#] /xms ) { 
        print "$input\n";
    }
}

