#!/usr/bin/env perl

# patch_pea_fa.pl -- Erich Schwarz <emsch@its.caltech.edu>, 1/13/2011.
# Purpose: fix .pea_fa files with uppercase residues (make all lowercase), and die loudly upon seeing 'n' or '.' residues.

use strict;
use warnings;

my $input = q{};

while (my $input = <>) { 
    chomp $input;
    if ( $input !~ /\A >/xms ) { 
        $input =~ tr/[A-Z]/[a-z]/;
        if ( $input =~ /[n]|[.]/xms ) { 
            die "PE-Assembler cannot accept 'n' or '.' residues!\n";
        }
    }
    print "$input\n";
}


