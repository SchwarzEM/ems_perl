#!/usr/bin/env perl

# fastq_trim_both.pl -- Ross Hall <rossh@unimelb.edu.au>, 7/1/2010.
# Purpose: "... reads sequences from a FASTQ file and outputs only first trim_length bases in each sequence.  Trims both sequence and quality."

# Usage: fastq_trim_both.pl input_FASTQ_file left_trim right_trim

use strict;
use warnings;

if ( @ARGV != 3 ) {
    print STDERR "fastq_trim_both.pl input_FASTQ_file left_trim right_trim \n";
    exit(1);
}

my $infile     = shift;
my $trim_left  = shift;
my $trim_right = shift;

open( INFILE, $infile ) || &ErrorMessage( "Cannot open file " . $infile );

my $line;
while (<INFILE>) {
    $line = $_;
    chomp($line);
    if ( $line =~ /^@/ && $line =~ /[0-9]/ ) {
        print "$line\n";
    }
    elsif ( $line =~ /^\+/ ) {
        print "$line\n";
    }
    else {
        my $sub1 = substr( $line, $trim_left );
        my $sub2 = substr( $sub1, 0, -$trim_right );
        print "$sub2\n";
    }
}

close(INFILE);

sub ErrorMessage {
    my $msg = shift;
    print "Fatal error: $msg\n";
    exit(1);
}

