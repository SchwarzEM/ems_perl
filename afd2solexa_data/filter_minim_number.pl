#!/usr/bin/env perl

# filter_minim_number.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/1/2009; minor edits 2/13/2010.
# Purpose: given a big space-delimited table with numbers in it, set minimum positive threshold value for *any* number seen.

use strict;
use warnings;

my $threshold = shift @ARGV;
if (! $threshold) { 
    die "Format: ./filter_minim_number.pl  [numerical threshold over 0]  [STDIN|file(s)]\n";
}
if (! ( $threshold > 0 ) ) { 
    die "Threshold $threshold was not a number above 0!\n";
}

while (my $input = <>) { 
    chomp $input;
    # The default decision is to print a line, but a single 'bad' value can reverse that.
    my $to_print = 1;

    # My own solution to the surprisingly tricky problem: *how* do you define a number?
    # Basically, I require that it be pure digits, or digits\.digits.
    # Note that I could also just stick in a subroutine requiring either >= or < threshold...
    my @data = grep { $_ !~ /[^\d|^\.]/xms } 
               grep { $_ =~ /\d\z/xms }
               grep { $_ =~ /\A\d/xms } 
               split /\s+/, $input;

    my $rpkm;
    foreach my $rpkm (@data) { 
        # Die loudly if a nonnumber did make it through.
        if ( (! ( $rpkm >= $threshold ) ) and (! ( $rpkm < $threshold ) ) ) { 
            die "Can't interpret datum $rpkm in input line $input\n";
        }
        # Finally: a single low value calls off the line printing.
        if ( $rpkm < $threshold ) { 
            $to_print = 0;
        }
    }
    if ($to_print) { 
        print "$input\n";
    }
}

