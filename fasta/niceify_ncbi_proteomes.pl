#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

# Sample input header:
# >gi|390331338|ref|XP_003723254.1| PREDICTED: uncharacterized protein LOC100887917 [Strongylocentrotus purpuratus]

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        if ( $input =~ /\A > (gi \| \d+ \| [a-z]+ \| ([^\|\s]+) \| [^\|\s]* \s+ .+ \[ ([A-Z]) [a-z]+ [ ] ([a-z]{3}) [a-z]* \] \s*) \z/xms ) { 
            my $full_header = $1;
            my $prot        = $2;
            my $sp_init     = $3;
            my $sp_lett3    = $4;
            my $prefix = $sp_init . $sp_lett3 . q{_};
            $prot = $prefix . $prot ;
            print ">$prot  $full_header\n";
        }
        else { 
            die "Can't parse header: $input\n";
        }
    }
    else { 
        print "$input\n";
    }
}

