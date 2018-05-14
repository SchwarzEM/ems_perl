#!/usr/bin/env perl

# bed2sizes.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/19/2010.
# Purpose: from a BED file, get an ascending table listing nt sizes (col. 1) and their frequency (col. 2).

use strict;
use warnings;

my $nt1   = q{};
my $nt2   = q{};
my $span1 = q{};

my %sizes = ();

while (my $input = <>) { 
    chomp $input;
    $nt1   = q{};
    $nt2   = q{};
    $span1 = q{};
    if ( $input =~ /\A \S+ \t (\d+) \t (\d+) \t \S+/xms ) { 
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
        $sizes{$span1} += 1;
    }
}

my @initial_spans = sort { $a <=> $b } keys %sizes;
my $top_span = $initial_spans[-1];

foreach my $span2 ( 1 .. $top_span ) { 
    if (exists $sizes{$span2}) { 
        print $sizes{$span2};
    }
    else { 
        print "0";
    }
    print "\t$span2\n";
}

