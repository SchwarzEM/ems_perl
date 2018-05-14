#!/usr/bin/env perl

# bed2size_slurry.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/19/2010.
# Purpose: from a BED file, get a single ascending column listing nt sizes, one by one (relying on Excel/oocalc/etc. to get a later frequency histogram).

use strict;
use warnings;

my $nt1   = q{};
my $nt2   = q{};
my $span1 = q{};

my @sizes = ();

while (my $input = <>) { 
    chomp $input;
    $nt1   = q{};
    $nt2   = q{};
    $span1 = q{};
    if (      ( $input =~ /\A \S+ \t (\d+) \t (\d+) \t \S+/xms ) 
          and ( $input !~ /\A track \s+ name= /xms                 ) ) { 
        $nt1 = $1;
        $nt2 = $2;
        if (! ( $nt1 < $nt2) ) { 
            die "Non-ascending coords. in: $input\n";
        }
        if ( (! ( $nt1 > 0 ) ) or (! ( $nt2 > 0 ) ) ) { 
            die "Non-positive coords. in: $input\n";
        } 
        # N.B.: this is assuming that the BED format doesn't actually give the first 
        #     nucleotide of the *element*, but of the last nucleotide before it starts!
        #     I suspect this is true because correcting by +1 nt gives me sizes that appear +1 nt too large.
        $span1 = $nt2 - $nt1;
        push @sizes, $span1;
    }
}

@sizes = sort { $a <=> $b } @sizes;

foreach my $span2 ( @sizes ) { 
    print "$span2\n";
}

