#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input ;
    if ( $input =~ /\A ((\S+ \/ ([^\s\/]+)) _1\.filt1\.fq\.gz) \t (\S*\/dbs\/(\S+)) \z/xms ) {
        my $input_1   = $1;
        my $full_stem = $2;
        my $file_stem = $3;
        my $db        = $4;
        my $db_stem   = $5;

        my $input_2  = $full_stem . '_2.filt1.fq.gz';

        $file_stem = "$file_stem.filt1";

        my $sam = "$file_stem.bwa-mem.$db_stem.sam";
        my $log = "$file_stem.bwa-mem.$db_stem.log.txt";
        my $err = "$file_stem.bwa-mem.$db_stem.err.txt";

        my $command = "bwa mem -v 3 -t 8 $db $input_1 $input_2 -o $sam 1>$log 2>$err ;";
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

