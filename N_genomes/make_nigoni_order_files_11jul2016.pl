#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $output = '/dev/null';

open my $OUTPUT, '>', $output;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A [>] (\S+) /xms ) { 
        my $bri_chr = $1;

        $output = "$bri_chr.order.txt";
        $output = safename($output);

        close $OUTPUT;
        open $OUTPUT, '>', $output;
    }
    elsif ( $output and ( $input =~ /\t ([+]|[-]) \t (\S+) \z/xms ) ) {
        my $ori = $1;
        my $seq = $2;
        print $OUTPUT "$ori$seq\n";
    }
    else {
        die "For output \"$output\", cannot parse input line: $input\n";
    }
}

close $OUTPUT;

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

