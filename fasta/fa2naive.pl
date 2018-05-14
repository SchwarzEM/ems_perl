#!/usr/bin/env perl

# fasta2naive.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/5/2009.
# Purpose: convert real FASTA files into naive one-line format (needed by RNAfold, et al.).

use strict;
use warnings;

my $reading = 0;

while (my $input = <>) { 
    if ( $input =~ /\S/xms ) { 
        if ( $input !~ /\A > \S/xms ) { 
            $input =~ s/\s+\z//;
            print "$input";
        }
        elsif ( $input =~ /\A > \S/xms ) {
            print "\n" if $reading;
            print $input;
            $reading = 1;
        }
        elsif ( $input =~ /\A > /xms ) { 
            chomp $input;
            die "Can't parse: $input\n";
        }
    }
}

# Last line, end of file:
print "\n";

