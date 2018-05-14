#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %type_count = ();

while (my $input = <>) {
    chomp $input;
    $type_count{$input}++;
}

my @inputs = sort { $type_count{$b} <=> $type_count{$a} }  keys %type_count;

foreach my $input (@inputs) {
    print "$input\t$type_count{$input}\n";
}


