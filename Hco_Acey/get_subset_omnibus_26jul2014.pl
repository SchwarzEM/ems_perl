#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $target_fasta = $ARGV[0];
my $full_headers = $ARGV[1];

my %seen = ();

open my $TARGET, '<', $target_fasta;
while (my $input = <$TARGET>) {
    chomp $input;
    if ( $input =~ /\A > ([^\s\/]+) \/ \d+ [-] \d+ /xms ) { 
        my $protein = $1;
        $seen{$protein} = 1;
    }
}
close $TARGET;

open my $HEADERS, '<', $full_headers;
while (my $input = <$HEADERS>) {
    chomp $input;
    if ( $input =~ /\A > (\S+) /xms ) {
        my $protein = $1;
        if ( exists $seen{$protein} ) {
            print "$input\n";
        }
    }
}
close $HEADERS;


