#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

# This script has one purpose: to enforce clean reformatting of read-1 FASTQ header lines from an unreliable source.  If a FASTQ header line fails to meet very basic criteria, the script dies, and thus saves me a lot of trouble later.

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A [@] /xms ) {
        if ( $input =~ /\A [@] \S+ [ ] \d: /xms ) {
            $input =~ s/\A([@]\S+)[ ]\d:/$1\#0\/1 1:/;
            print "$input\n";
        }
        else {
            die "Misformatted FASTQ header line: $input\n";
        }
    }
    else {
        print "$input\n";
    }
}

