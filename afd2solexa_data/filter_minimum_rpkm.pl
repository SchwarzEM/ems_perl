#!/usr/bin/env perl

# filter_minimum_rpkm.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/13/2010.
# Purpose: weed out RPKM values (in third column) either below minimum threshold, or not above greater-than threshold.

use strict;
use warnings;
use Getopt::Long;

my $minimum;
my $greater_than;

GetOptions ( 'minimum:f'      => \$minimum,
             'greater_than:f' => \$greater_than, );

if (    ( (! defined $minimum ) and (! defined $greater_than) ) 
     or ( ( defined $minimum ) and ( defined $greater_than ) ) ) { 
    die "Format: filter_minimum_rpkm.pl",
        " --minimum|-m [minimum RPKM allowed] <or>",
        " --greater_than|-g [maximum RPKM *not* allowed]",
        " <input stream/files>\n",
        ; 
}

if ( defined $minimum )  { 
    if (! ( $minimum > 0 ) ) { 
        die "Minimum allowed RPKM value $minimum was not a number above 0.00!\n";
    }
}

if ( defined $greater_than ) { 
    if (! ( $greater_than >= 0 ) ) {
        die "Maximum disallowed RPKM value $greater_than was not a number equal to or above 0.00!\n";
    }
}

while (my $input = <>) { 
    chomp $input;
    my $rpkm;
    my $threshold;
    if ( ( $input =~ /\A \S+ \s+ \S+ \s+ (\S+) /xms ) and ( $input !~ /\A \# /xms ) ) { 
        $rpkm = $1;
        if ( ( defined $minimum ) and ( $rpkm >= $minimum ) ) { 
            print "$input\n";
            $threshold = $minimum;
        }
        elsif ( ( defined $greater_than ) and ( $rpkm > $greater_than ) ) { 
            print "$input\n";
            $threshold = $greater_than;
        }
    }
    # Sanity-check $rpkm to make sure it's a number.
    if (     ( defined $rpkm            ) 
         and ( defined $threshold       ) 
         and (! ( $rpkm >= $threshold ) ) 
         and (! ( $rpkm < $threshold  ) ) ) { 
        die "Can't parse input line: $input\n";
    }
}

