#!/usr/bin/env perl

# coordfilt_ali_sites.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/2/2008.
# Purpose: extract only Ali sites within a defined range.

use strict;
use warnings;
use Getopt::Long;

my $coords = q{};
my $chromo = q{};
my $nt1    = q{};
my $nt2    = q{};

GetOptions ( "coords=s" => \$coords );

if ( (! $coords ) 
      or ( $coords !~ / \A (I|II|III|IV|V|X) : \d+ \.\. \d+ \z /xms) ) { 
    die "Format: ./coordfilt_ali_sites.pl",
        "  [-c /--coords=]CHR:X..Y  [input file]\n",
        ;
}
if ( $coords =~ / \A (I|II|III|IV|V|X) : (\d+) \.\. (\d+) \z /xms) { 
    $chromo = $1;
    $nt1    = $2;
    $nt2    = $3;
}

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\d+) \t (\d+) /xms) {
        my $curr_chr = $1;
        my $curr_nt1 = $2;
        my $curr_nt2 = $3;
        if ( ( $curr_chr eq $chromo  ) 
             and ( $nt1 <= $curr_nt1 ) 
             and ( $curr_nt2 <= $nt2 ) ) {
            print "$input\n";
        }
    }
}

