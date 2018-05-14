#!/usr/bin/env perl

# simple_fastq2fa.pl -- Erich Schwarz <ems394@cornell.edu>, 1/20/2016.
# Purpose: very simple-minded FastQ to FastA converter; does not deal with non-stereotypically formatted FastQ but should not jam on third-line variants.

use strict;
use warnings;
use autodie;

my $i = 0;

while (my $input = <>) {
    chomp $input;
    $i++;
    if ( ($i % 4) == 1 ) {
        if ( $input =~ /\A [@] /xms ) { 
            $input =~ s/\A[@]/>/;
            print "$input\n"; 
        }
        else {
            die "Cannot parse putative first-line header: $input\n";
        }
    }
    elsif ( ($i % 4) == 2 ) {
        print "$input\n";
    }
}

