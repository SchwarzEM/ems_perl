#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = "\tcondition";

while (my $input = <>) {
    chomp $input;
    my @data_types = split "\t", $input;
    foreach my $data (@data_types) { 
        if ( $data =~ /\A (\S+) _rep\d+_reads \z/xms ) {
            my $condition = $1;
            print "$header\n" if $header;
            $header = q{};
            print "$data\t$condition\n";
        }
    }
}

