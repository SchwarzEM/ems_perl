#!/usr/bin/env perl

use strict;
use warnings;

my $val       = q{};
my $consensus = q{};

while (my $input = <>) { 
     chomp $input;
     if ( $input =~ / \A (?: \d+ (?: \d*\.\d+){0,1} \s+){5} (?: \S+ \s+){2} (\S+) /xms ) { 
         $val = $1;
         $consensus .= $val;
     }
}
print "\nConsensus of $ARGV: $consensus\n";



