#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $i = 4;
my $j = 0;

while (my $input = <>) {
    chomp $input;
    $i++;
    $j = int ($i / 5);
    $j = sprintf "%02i", $j;

    my $outdir = 'rhabditella_' . $j;

    if ( int ($i / 5) == ($i / 5) ) {
        print "\n    mkdir $outdir ;\n";
    }
    print "    mv -i rhabditella/$input $outdir ;\n";
}

