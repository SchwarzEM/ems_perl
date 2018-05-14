#!/usr/bin/perl

# count_integr8_ids.pl, Erich Schwarz <emsch@its.caltech.edu>, 3/2/06.

use strict;
use warnings;
use File::Basename;

my %count = ();

while (<>) { 
    chomp (my $input = $_);
    if ($input =~ /^ID\s+/) { 
        $count{basename($ARGV)}++;
    }
}
    

foreach my $f (sort { $count{$a} <=> $count{$b} } keys %count) { 
    print "$count{$f}\t$f\n";
}
