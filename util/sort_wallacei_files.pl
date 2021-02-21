#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $i = 3;
my $j = 0;

while (my $input = <>) {
    chomp $input;
    $i++;
    $j = int ($i / 4);
    $j = sprintf "%02i", $j;

    my $outdir = 'wallacei_' . $j;

    if ( int ($i / 4) == ($i / 4) ) {
        print "\n    mkdir $outdir ;\n";
    }
    print "    mv -i wallacei/$input $outdir ;\n";
}

