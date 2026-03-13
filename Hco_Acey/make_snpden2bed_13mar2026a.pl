#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while ( my $infile = <> ) {
    chomp $infile;
    if ( $infile =~ /\A (\S+) \. snpden \z/xms ) {
        my $stem = $1;
        my $outfile = "$stem.bed";
        my $err     = "$stem.err";
        print '$PROJECT/ems_perl/Hco_Acey/snpden2bed_13mar2026.pl ';
        print "$infile 1>$outfile 2>$err ;\n";
    }
    else {
        die "Cannot parse infile: $infile\n";
    }
}

