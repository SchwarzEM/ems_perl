#!/usr/bin/env perl

# sto2rna.pl -- Erich Schwarz <ems394@cornell.edu>, 7/21/2013.
# Purpose: given an ACGT *.sto alignment to be used for INFERNAL, produce the same alignment, but with ACGU.  May or may not be actually needed, but good for peace of mind.

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( ( $input !~ /\A[#]/xms ) and ( $input !~ /\A\//xms ) and ( $input =~ /\A \S /xms ) ) { 
        if ( $input =~ /\A (\S+ \s+) (\S+) \z/xms ) {
            my $header   = $1;
            my $sequence = $2;
            $sequence =~ s/T/U/g;
            $sequence =~ s/t/u/g;
            print "$header$sequence\n";
        }
        else { 
            die "Can't parse input: $input\n";
        }
    }
    else { 
        print "$input\n";
    }
}

