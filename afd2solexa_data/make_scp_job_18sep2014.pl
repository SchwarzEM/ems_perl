#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = '#!/bin/bash' . "\n\n";
my @lines  = ();
my $footer = "\n";


while (my $input = <>) {
    chomp $input;
    my $command = "    scp -rp $input " . 'emsch@hpcc.msu.edu:/mnt/home/emsch/work/AFD_etc/rnaseq_data/rnaseq_reads/modENCODE ;' . "\n";
    push @lines, $command;
}

if (@lines) {
    print $header;
    print @lines;
    print $footer;
}    

