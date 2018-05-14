#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;

    $input =~ s/_rep(\d)/ rep. $1/;
    $input =~ s/atml1_(\d)/atml1-$1/;
    $input =~ s/__/::/;
    $input =~ s/lgo_(\d)/lgo-$1/;
    $input =~ s/LGO_/LGO /;
    $input =~ s/FLAG_/FLAG:/;
    $input =~ s/ATML1_/ATML1 /;

    print "$input\n";
}


