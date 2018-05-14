#!/usr/bin/env perl

# extract_fastq_namelist.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/2/2010.
# Purpose: quickly extract a simple list of names from a FASTQ file: ignore redundant names, but do reject obviously misformatted files.

use strict;
use warnings;

my $line_no = 0;
my $i       = 0;
my %names   = ();

while (my $input = <>) { 
    chomp $input;
    $i++;
    $line_no = ($i % 4);
    if ( $line_no == 1 ) { 
        if ( $input =~ /\A @ (\S+) /xms ) { 
            my $seq = $1;
            print "$seq\n";
        }
        else { 
            die "Can't parse putative header line: $input\n";
        }
    }
}

