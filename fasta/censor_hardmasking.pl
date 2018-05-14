#!/usr/bin/env perl

# censor_hardmasking.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/20/2008.
# Purpose: filter/censor 'n', 'N', 'x', or 'X' residues from DNA/RNA FASTAs.

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > \S+ /xms ) { 
        print "$input\n";
    }
    else { 
        if ( $input =~ / [^acgntuxACGNTUX] /xms ) { 
            die "Detected non-DNA/RNA-like residues in: $input\n";
        }
        $input =~ s/[nxNX]//g;
        # If the *entire* line's gone, don't bother printing:
        if ( $input =~ /\S/xms ) { 
            print "$input\n";
        }
    }
}

