#!/usr/bin/env perl

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A>/xms ) { 
        if ( $input =~ /\A > (\S+) /xms ) { 
            my $seq  = $1;
            my $gene = $seq;
            $gene =~ s/[a-z]\z//;
            if ( $input =~ / locus: (\S+) /xms ) { 
                $gene = $1;
            }
            $gene = 'Cel_' . $gene;
            print ">$gene\n";
        }
        else {
            die "Can't parse header: $input\n";
       }
    }
    else { 
        print "$input\n";
    }
}

