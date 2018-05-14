#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A \S+ \t \S+ \t (\d+) \t (\d+) \z/xms ) { 
        my $start_nt = $1;
        my $end_nt   = $2;
        $data_ref->{'start_nt'}->{$start_nt}->{'end_nt'}->{$end_nt}->{'match'}->{$input} = 1;
    }
    else {
        die "Cannot format input line: $input\n";
    }
}

my @starts = sort { $a <=> $b } keys %{ $data_ref->{'start_nt'} };
foreach my $start_nt (@starts) {
    my @ends = sort { $a <=> $b } keys %{ $data_ref->{'start_nt'}->{$start_nt}->{'end_nt'} };
    foreach my $end_nt (@ends) {
        my @matches = sort keys %{ $data_ref->{'start_nt'}->{$start_nt}->{'end_nt'}->{$end_nt}->{'match'} };
        foreach my $match (@matches) {
            print "$match\n";
        }
    }
}


