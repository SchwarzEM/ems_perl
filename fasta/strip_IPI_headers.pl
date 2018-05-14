#!/usr/bin/env perl

# strip_IPI_headers.pl -- Erich Schwarz <emsch@its.caltech.edu>, 7/10/2009.
# Purpose: reformat '>IPI:IPI00000001.2 ...' to just '>IPI00000001.2'; die loudly if nonunique name arises. 

use strict;
use warnings;

my %seen = ();

while (my $input = <>) { 
    if ( $input !~ /\A > /xms ) { 
        print $input;
    }
    elsif ( $input =~ /\A > IPI: ( [^(?:\s|\|)]+ ) \| /xms ) { 
        my $protein = $1;
        if ($seen{$protein}) { 
            die "Protein name $protein nonredundant!\n";
        }
        print ">$protein\n";
        $seen{$protein} = 1;
    }
    elsif ( $input =~ /\A > /xms ) { 
        chomp $input;
        die "Can't parse header: $input\n";
    }
}

