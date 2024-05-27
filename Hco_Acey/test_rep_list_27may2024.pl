#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

while (my $input = <>) {
    chomp $input;
    $input =~ s/\\//g;
    $input =~ s/\A\s+//;
    $input =~ s/\s+\z//;
    my @list = split /\s+/, $input;
    @list = uniq(@list);
    print "@list\n";
}
