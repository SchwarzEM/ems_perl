#!/usr/bin/env perl

use strict;
use warnings;

while (my $input = <>) {
    chomp $input;
    my $new_name = $input . '.orig_revcomp';
    if (-e $new_name) {
        die "Try a different renaming for $input than $new_name\n";
    }
    print "    mv -i $input $new_name;\n";
}


