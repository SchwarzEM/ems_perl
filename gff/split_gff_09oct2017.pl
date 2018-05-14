#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( ( $input !~ /\A [#] /xms ) and ( $input =~ /\A (\S+) /xms ) ) { 
        my $seq = $1;
        my $gff = "$seq.gff";
        open my $GFF, '>>', $gff;
        print $GFF "$input\n";
        close $GFF;
    }
}

