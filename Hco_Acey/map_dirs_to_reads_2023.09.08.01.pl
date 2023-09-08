#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while ( my $input = <> ) {
    chomp $input;
    if ( $input =~ /\A \/ (N\d+) \/ (\S+) \. \z/xms ) {
        my $dir   = $1;
        my $stem  = $2;
        my $file1 = "$stem.fastq.gz";
        my $file2 = $file1;

        $file2 =~ s/_R1_/_R2_/g;
        print "$file1\t$dir\n";
       	print "$file2\t$dir\n";
    }
    else {
        warn "Cannot parse: $input\n";
    }
}
