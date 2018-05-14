#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;

my @input_files = ();

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) /xms ) {
        my $file = $1;
        if (! -r $file) {
            die "Cannot find readable copy of this file: $file\n";
        }
        if ( ( $file !~ /FLAG_ATML1/xms ) and ( $file !~ /atml1_4/xms ) ) { 
            push @input_files, $file;
        }
    }
}

my @output_files = sort { basename($a) cmp basename($b) } @input_files;

foreach my $output_file (@output_files) {
    print "$output_file\n";
}

