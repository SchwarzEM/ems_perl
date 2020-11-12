#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input ;
    if ( $input =~ /\A ((\S+ \/ ([^\s\/]+)) _1\.filt1\.fq\.gz) \t (\S*\/dbs\/(\S+)\.hisat2) \z/xms ) {
        my $input_1   = $1;
        my $full_stem = $2;
        my $file_stem = $3;
        my $db        = $4;
        my $db_stem   = $5;

        my $input_2  = $full_stem . '_2.filt1.fq.gz';

        $file_stem = "$file_stem.filt1";

        my $sam = "$file_stem.hisat2.$db_stem.sam";
        my $log = "$file_stem.hisat2.$db_stem.log.txt";
        my $err = "$file_stem.hisat2.$db_stem.err.txt";

        my $command = "hisat2 -p 8 -x $db -1 $input_1 -2 $input_2 -S $sam 1>$log 2>$err ;";
        print "$command\n";

        if (! -r $input_1 ) {
            die "Can't read input 1: $input_1\n";
        }
        if (! -r $input_2 ) {
            die "Can't read input 2: $input_2\n";
        }

    }
    else {
        die "Cannot parse: $input\n";
    }
}

