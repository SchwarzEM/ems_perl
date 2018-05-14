#!/usr/bin/env perl

# compress_script_copy.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/24/2010.
# Purpose: given human-readable script lines in electronic ASCII text, rework into almost-runnable script text.

use strict;
use warnings;

my $script_text;

while ( my $input = <> ) { 
    chomp $input;
    $script_text .= $input;
}

$script_text =~ s/\\//g;
$script_text =~ s/[ ]+/ /g;
$script_text =~ s/;/;\n/g;
$script_text =~ s{(\#!\/\S+perl)}{$1\n\n};
$script_text =~ s{(\#!\/\S+bash)}{$1\n\n};
$script_text =~ s/\n[ ]/\n    /g;
print $script_text;
print "\n";

