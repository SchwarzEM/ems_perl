#!/usr/bin/env perl

# reject_sheep_Hco_contigs_27aug2011.pl -- Erich Schwarz <emsch@caltech.edu>, 8/27/2011.
# Purpose: select which Hco contigs to remove either by >=50% match to sheep/cow or >=100-nt/>=5% match.

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

my $match   = q{};
my $percent = q{};

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A \S+ \t \S+ \t (\S+) \t (\S+) /xms ) { 
        $match   = $1;
        $percent = $2;
        $match =~ s/[,]//g;
        if ( looks_like_number($match) and looks_like_number($percent) ) { 
            if (    ( $percent >= 50                          ) 
                 or ( ( $match >= 100 ) and ( $percent >= 5 ) ) ) {  
                print "$input\n";
            }
        }
    }
}

