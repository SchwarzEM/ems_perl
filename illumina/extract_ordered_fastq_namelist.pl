#!/usr/bin/env perl

# extract_ordered_fastq_namelist.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/2/2010.
# Purpose: get an ordered list of names from a FASTQ file: violently reject redundant names or obviously misformatted files.

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
        if ( $input !~ /\A @ \S+ /xms ) { 
            die "Can't parse putative header line: $input\n";
        }
        if ( $input =~ /\A @ (\S+) /xms ) { 
            my $seq = $1;
            if ( exists $names{$seq} ) {
                die "Non-unique sequence name $seq in putative header line: $input\n";
            }
            $names{$seq} = 1;
        }
    }
}

foreach my $name (sort keys %names) { 
    print "$name\n";
}

