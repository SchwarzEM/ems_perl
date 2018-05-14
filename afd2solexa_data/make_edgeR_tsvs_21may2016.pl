#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $transform = q{| perl -ne ' s/\"\"/Gene/; s/[,]/\t/g; s/["]//g; print; ' >};

while (my $infile = <>) {
    chomp $infile;
    my $outfile = $infile;
    $outfile =~ s/\.csv\z/.tsv.txt/;
    if (-e $outfile) {
        die "Will not target an existing file for output: $outfile\n";
    }
    print "   cat $infile $transform $outfile ;\n";
}

