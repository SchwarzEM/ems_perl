#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A\/\S+\/(\S+)\.R1\.filt1\.fq\z/xms ) {
        my $prefix = $1;
        print "$prefix\t$input\n";
    }
}
