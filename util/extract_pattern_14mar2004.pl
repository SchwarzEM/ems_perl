#!/usr/bin/perl

# extract_pattern_14mar2004.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/14/04
# Purpose: extract file-defined pattern from text -- ./extract_pattern.pl  pattern_file  input_file  [output_file]

# Example 'pattern' file (must be one line and must have parentheses):
# WB:(\[.+\])\t

use strict;

my $pattern_file      = "";
my $input_file        = "";
my $output_file       = "";
my $pattern_to_match  = "";

die "Usage:  ./extract_pattern.pl  pattern_file  input_file  [output_file]\n" if (!$ARGV[1]);
die "Usage:  ./extract_pattern.pl  pattern_file  input_file  [output_file]\n" if ($ARGV[3]);

chomp ($pattern_file  = $ARGV[0]);
chomp ($input_file    = $ARGV[1]);

if ($ARGV[2]) {
    chomp($output_file  = $ARGV[2]);
}
elsif (!$ARGV[2]) {
    $output_file        = $input_file . ".output";
}

open (PATTERN, "$pattern_file")  || die "Can't open $pattern_file $!\n";
open (INPUT,   "$input_file")    || die "Can't open $input_file $!\n";
open (OUTPUT,  ">$output_file")  || die "Can't open $output_file $!\n";

chomp($pattern_to_match = <PATTERN>);

while (<INPUT>) {
    chomp(my $input_line = $_);
    if ($input_line =~ /$pattern_to_match/) {
        print OUTPUT "$1\n";
        # N.b.: $pattern_to_match MUST have .*\(.*\).*
    }
}

close PATTERN;
close INPUT;
close OUTPUT;
