#!/usr/bin/env perl

# regex_ify.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/22/2010.
# Purpose: convert plain text into its regex equivalent, as nearly as possible.

use strict;
use warnings;

my %equiv = ( q{ } => q{[ ]},
              q{.} => q{\.},
              q{-} => q{\-}, 
              q{|} => q{\|}, );

while (my $input = <>) { 
    chomp $input;
    my @in_chars  = split //, $input;
    my @out_chars = ();
    foreach my $char (@in_chars) { 
        if ( exists $equiv{$char} ) { 
            $char = $equiv{$char};
        }
        push @out_chars, $char;
    }
    my $output = join q{}, @out_chars;
    print "$output\n";
}

