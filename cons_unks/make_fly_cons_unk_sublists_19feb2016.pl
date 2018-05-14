#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $single_out = 'Drosophila_unique_cons_unks_Oct2014.tsv.txt';
my $mult_out   = 'Drosophila_multiple_cons_unks_Oct2014.tsv.txt';

$single_out = safename($single_out);
$mult_out   = safename($mult_out);

open my $SINGLE, '>', $single_out;
open my $MULT,   '>', $mult_out;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ / drosophila .+ drosophila /xms ) {
        print $MULT "$input\n";
    }
    elsif ( $input =~ / drosophila /xms ) {
        print $SINGLE "$input\n";
    }
}

close $SINGLE;
close $MULT;

sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}


