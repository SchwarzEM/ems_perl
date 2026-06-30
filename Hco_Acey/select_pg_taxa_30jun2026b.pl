#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while ( my $input = <> ) { 
    chomp $input;
    if ( $input =~ /\A Gene \t/xms ) {
        print "$input\n";
    }
    elsif ( $input =~ /\A [^\t]+ \t [^\t]+ \t ([^\t]+) \t ([^\t]*) \t ([^\t]*) /xms ) {
        my $pass_taxa   = $1;
        my $array_taxa  = $2;
        my $nonsyn_taxa = $3;
        if (     ( $pass_taxa =~ /Ilik2/xms ) 
             and ( $pass_taxa !~ /Aroian/xms ) 
             and ( $pass_taxa =~ /(Oita|obscurus)/xms ) 
             and ( $array_taxa eq 'Baylor; Keiser' )
           ) {
            print "$input\n";
        }
    }
    else {
        die "Cannot parse: $input\n";
    }
}

