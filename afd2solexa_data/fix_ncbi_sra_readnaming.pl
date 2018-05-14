#!/usr/bin/env perl

# fix_ncbi_sra_readnaming.pl -- Erich Schwarz <emsch@caltech.edu>, 3/2/2012.
# Purpose: given an NCBI paired-end FASTQ file with idiotic naming system, hack the names into something which I can reasonably parse with my fastq scripts.

use strict;
use warnings;

my $stem = q{};

my $head_pt1 = q{};
my $head_pt2 = q{};
my $head_pt3 = q{};

my $i = 0;
my $j = 0;

while (my $input = <>) { 
    chomp $input;
    $i++;
    $j = ( $i % 8 );

    if ($j == 1) {
        if ( $input =~ /\A [@] (SRR\d+\.\d+) [ ] (\S+) [ ] (length=\d+) \s* \z/xms ) {
            $head_pt1 = $1;
            $head_pt2 = $2;
            $head_pt3 = $3;
            $stem = $head_pt1 . '_' . $head_pt2;
            print '@', $stem , '#0/1', q{ } , "$head_pt3\n" ;
        } 
        else { 
            die "Can't parse header line: $input\n";
        }
    }

    if ( ($j == 2) or ($j == 4) or ($j == 6) or ($j == 0) ) { 
        print "$input\n";
    }

    if ($j == 5) {
        if ( $input =~ /\A [@] $head_pt1 [ ] $head_pt2  [ ] $head_pt3 \s* \z/xms ) {
            print '@', $stem , '#0/2', q{ } , "$head_pt3\n" ;
        }
        else {
            die "Can't parse header line: $input\n";
        }
    }

    if ( ($j == 3) or ($j == 7) ) {
        if ( $input =~ /\A [+] $head_pt1 [ ] $head_pt2  [ ] $head_pt3 \s* \z/xms ) {
            print '+', $stem ;
            print '#0/1' if ($j == 3);
            print '#0/2' if ($j == 7);
            print q{ } , "$head_pt3\n" ;
        }
        else {
            die "Can't parse header line: $input\n";
        }
    }
}

