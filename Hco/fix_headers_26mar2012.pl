#!/usr/bin/env perl

# fix_headers_26mar2012.pl -- Erich Schwarz <emsch@caltech.edu>, 3/26/2012.
# Purpose: rework erroneous headers of a genome assembly from 'aa' to 'nt' counts.

use strict;
use warnings;

my $front = q{};
my $back  = q{};

while (my $input = <>) { 
    chomp $input;
    if ( $input !~ /\A > /xms ) { 
        print "$input\n";
    }
    else { 
        if ( $input =~ /\A > (\S+ \s+ \S+ [ ]) aa ([ ] .*) \z/xms ) { 
            $front = $1;
            $back  = $2;
            print '>', $front, 'nt', $back, "\n";
         }
        else { 
            die "Can't parse header line: $input\n";
        }
    }
}

