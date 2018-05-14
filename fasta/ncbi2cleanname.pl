#!/usr/bin/env perl

# ncbi2cleanname.pl -- Erich Schwarz <emsch@caltech.edu>, 11/30/2011.
# Purpose: given something like est_others from NCBI with names that choke other programs reading it (e.g., getorf/EMBOSS), convert its names into unique derivatives of the foo|XXX|zork|YYY naming NCBI uses, while keeping the original name and other text in the header.

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > ( (?: [^\s\|]+ \|)+ ([^\s\|]+) \| \s .* )\z/xms ) {  
        my $header = $1;
        my $name   = $2;
        $name = safename($name);
        print ">$name  $header\n";
    }
    else { 
        print "$input\n";
    }
}

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

