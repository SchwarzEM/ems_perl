#!/usr/bin/env perl

use strict;
use warnings;

# PF10573|UPF0561	11.5009	4 genes	c_elegans|WBGene00010807|M01E5.4; hu

my $i = 1;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \t \S+ \t [^\t]* \t ([^\t]+) \z/xms ) { 
        my $front = $1;
        my $genes = $2;
        $genes =~ s/\A.*(c_elegans.+?);.*/$1/;
        $genes =~ s/c_elegans\|//g;
        print "$i\t$front\t$genes\n";
        $i++;
    }
    else { 
        die "Can't parse input: $input\n";
    }
}


