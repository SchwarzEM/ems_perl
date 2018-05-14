#!/usr/bin/env perl

# scnd2last_ali_motlist.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/2/2008.
# Purpose: unique genome sites w/ >= 1 sources/, printed in genomic order.

use strict;
use warnings;

my $coords = q{};
my $source = q{};
my $rawseq = q{};
my $coords2nos_ref;
my $coords2data_ref;
my %coords2rawseq = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input !~ / \A \S+ \t \d+ \t \d+ \t \S+ \t [ACGT]+ \z /xms ) { 
        warn "Could not parse input line: $input\n";
    }
    if ( $input =~ / \A (\S+ \t \d+ \t \d+) 
                        \t (\S+) 
                        \t ([ACGT]+) 
                     \z /xms ) { 
        $coords = $1;
        $source = $2;
        $rawseq = $3;
        my @coord_list = split "\t", $coords;
        $coords2rawseq{$coords} = $rawseq; 
        $coords2data_ref->{$coords}->{$source} = 1;
        $coords2nos_ref->{$coords} = \@coord_list;
    }
}

my @sorted_coords = 
    sort {     $coords2nos_ref->{$a}->[1] 
           <=> $coords2nos_ref->{$b}->[1] }
    sort {     $coords2nos_ref->{$b}->[2]
           <=> $coords2nos_ref->{$a}->[2] }
    keys %{ $coords2nos_ref };

foreach my $sort_coord (@sorted_coords) { 
    my @sources = sort keys %{ $coords2data_ref->{$sort_coord} };
    my $sourceline = join '|', @sources;
    print $sort_coord,
          "\t",
          $sourceline,
          "\t",
          $coords2rawseq{$sort_coord},
          "\n",
          ;
}

