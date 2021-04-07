#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %names = ();

while (my $input = <>) {
    chomp $input;
    $names{$input}++;
    my $suffix = $names{$input};
    print "$input$suffix\n";
}


