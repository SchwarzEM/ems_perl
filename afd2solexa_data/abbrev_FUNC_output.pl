#!/usr/bin/env perl

# abbrev_FUNC_output.pl -- Erich Schwarz <emsch@its.caltech.edu>, 7/23/2011.
# Purpose: get rid of most FUNC output columns, which end readers almost certainly do not need.

use strict;
use warnings;

my $front = q{};
my $back  = q{};

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A ( [^\t]+ \t [^\t]+ \t [^\t]+ ) \t .* (\t [^\t]* \S) \s* \z/xms ) { 
        $front = $1;
        $back  = $2;
        $input = $front . $back;
        print "$input\n";
    }
    else { 
        die "Can't parse: $input\n";
    }
}

