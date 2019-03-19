#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

# Sample input:
# #=GF ID   1-cysPrx_C
# #=GF AC   PF10417.9

my $acc  = q{};
my $name = q{};

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A [#]=GF [ ] ID [ ]+ (\S+)/xms ) {
        $name = $1;
    }
    elsif ( $name and ( $input	=~ /\A [#] =GF [ ] AC [ ]+ (\S+) /xms ) ) {
         $acc = $1;
         if ( $acc =~ /\A (PF\d+) \.\d+ \z/xms ) {
       	     $acc = $1;

             print "$acc\tpfam\|$acc\|$name\n";

             $acc  = q{};
             $name = q{};
       	 } 
         else {
             die "Cannot parse Pfam acc $acc in: $input\n";
         }
    }
}
