#!/usr/bin/env perl

# fastq2fa_simple.pl -- Erich Schwarz <emsch@caltech.edu>, 3/3/2012.
# Purpose: given a very simple input of FASTQ -- four lines per stanza -- turn the first two lines into FASTA text, and discard the remaining two.

use strict;
use warnings;

my $i      = 0;
my $j      = 0;
my $header = q{};

while (my $input = <>) { 
    chomp $input;
    $i++;
    $j = ( $i % 4 );
    if ( $j == 1 ) { 
        if ( $input =~ /\A [@] (\S+ .*) \z/xms ) {
            $header = $1;
            print '>' . $header . "\n";
        }
        else { 
            die "Can't parse putative header line: $input\n";
        }
    }
    if ( $j == 2 ) { 
        if ( $input =~ /\A [ACGTNacgtn]+ \z/xms ) { 
            print "$input\n";
        }
        else {
            die "Can't parse putative DNA sequence line: $input\n";
        }
    }
}

