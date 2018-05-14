#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $new_header = 'batch';

my %num2word = (
    '1' => 'one',
    '2' => 'two',
    '3' => 'three',
);

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A\t/xms ) { 
        $input =~ s/\A\t/\t$new_header\t/;
    }
    elsif ( $input =~ /\A (\S+_rep(\d+)_reads)\t(.*)\z/xms ) { 
        my $head = $1;
        my $rep  = $2;
        my $tail = $3;
        $rep     = $num2word{$rep};
        $input = "$head\t$rep\t$tail";
    }
    print "$input\n";
}

