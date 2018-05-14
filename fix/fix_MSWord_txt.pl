#!/usr/bin/perl

# fix_MS_returns.pl -- Erich Schwarz <emsch@its.caltech.edu>, 4/1/2009
# Purpose: clean up "\r" in defective Word-to-text, either inline or stream.

use strict;
use warnings;
use Getopt::Long;

my $inline = q{};

GetOptions ( 'inline|backup' => \$inline, ); 

if ( $inline ) { 
    $^I = '.bak';
}

while (my $input = <>) { 
    $input =~ s/\r/\n/g;
    print "$input";
}

