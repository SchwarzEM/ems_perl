#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

print '#!/bin/bash', "\n\n";

while (my $input = <>) { 
    chomp $input;
    $input =~ s/\A\s+//;
    $input =~ s/\s+\z//;
    my $basename = basename $input;
    if ( $basename =~ /\A (\S+) \.all\.iprscan\.raw\.txt \z/xms ) { 
        my $stem = $1;
        my $output = $stem . '_total_iprscan_raw.txt';
        $output    = safename($output);
        print "    cut -f 4-6,9 $input > $output;\n";
    }
    else { 
        die "Can't parse $input\n";
    }
}

print "\n";

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


