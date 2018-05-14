#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    $input =~ s#/mnt/ls15/scratch/users/emsch/work_rsync/2015/caenogens/nigoni/pbjelly#/mnt/ls15/scratch/users/emsch/work_rsync/2015/caenogens/tropicalis/pbjelly_25x#g;
    $input =~ s/nigoni/tropicalis/g;
    $input =~ s/Nigoni/Tropicalis_25x/g;
    $input =~ s/2015\.11\.11/25x_2016.07.29/g;
    print "$input\n";
}

