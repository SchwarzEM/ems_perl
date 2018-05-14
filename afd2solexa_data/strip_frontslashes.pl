#!/usr/bin/env perl

# strip_frontslashes.pl -- Erich Schwarz <emsch@its.caltec.edu>, 11/23/2010.
# Purpose: remove dorky ACeDB-style frontslashes from annotation text.  Archived as a script when perl -ne proved unable to grok the commands in this script!

use strict;
use warnings;

while (my $input = <>) { 
    $input =~ s/\\//g;
    print $input;
}

