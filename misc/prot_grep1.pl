#!/usr/bin/perl

# prot_grep1.pl, Erich Schwarz <emsch@its.caltech.edu>, 3/3/06.
# Purpose: make substitute for apparently buggy GNU sort.

use strict;
use warnings;

my %proteome_count = ();

while (<>) { 
    chomp(my $input = $_);
    unless ($input =~ /^#\s/) { 
        if ($input =~ /(.+)\t(\d+)$/) { 
            $proteome_count{$1} = $2;
        }
    }
}

foreach my $p (sort { $proteome_count{$b} 
                      <=> $proteome_count{$a} } 
                    keys %proteome_count) { 
     print "$proteome_count{$p}\t$p\n";
}

