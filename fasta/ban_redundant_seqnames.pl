#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $print = 1;
my %seen  = ();

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A > /xms ) {
        if ( $input =~ /\A > (\S+) /xms ) { 
            my $seqname = $1;
            if (! $seen{$seqname} ) {
                $print = 1;
                print "$input\n";
                $seen{$seqname} = 1;
            }
            else { 
                $print = 0;
            }
        }
        else { 
            die "Can't parse empty header: $input\n";
        }
    }
    else { 
        print "$input\n" if $print;
    }
}
