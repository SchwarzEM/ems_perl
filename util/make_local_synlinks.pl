#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

while (my $input = <>) { 
    chomp $input;
    my $target = basename $input;
    print "    ln -s $input $target;\n";
}

