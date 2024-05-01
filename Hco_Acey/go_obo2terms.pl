#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $id   = q{};
my $desc = q{};

while (my $input = <> ) {
    # id: GO:0000001
    # name: mitochondrion inheritance
    chomp $input;
    if ( $input =~ /\A id: \s (GO:\d+) /xms ) {
        $id = $1;
    }
    elsif ( $input =~ /\A name: \s (.*\S) /xms ) { 
        $desc = $1;
        $desc =~ s/\A\s+//;
        print "$id\t$desc\n";
        $id   = q{};
        $desc = q{};
    }
}

