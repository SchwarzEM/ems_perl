#!/usr/bin/perl

# grab_min_nonzero_evalue.pl -- Erich Schwarz <emsch@its.caltech.edu>, 6/26/2006.
# Purpose: list non-zero E values (lowest up) from all-on-all BlastP, e.g., to properly calibrate OrthoMCL.

use strict;
use warnings; 

# sample input:
# Sequences producing significant alignments:                      (bits) Value
# 
# CBP04552                                                              670   0.0
# C40C9.1 CE24839 WBGene00006673 locus:twk-20 potassium channel pr...   262   5e-70

my $input = "";
my $take_e_value = 0;
my $e_value = 0;
my %e_value_list = ();

while (<>) { 
    chomp($input = $_);
    if (($take_e_value) and ($input =~ /^.+\d+\s+(\S+)\s*$/)) {
        $e_value = $1;
        if ($e_value =~ /^\d+e-\d+$/) { 
            $e_value_list{$e_value} = 1;
        }
        elsif ($e_value =~ /^e-\d+$/) {
            $e_value = "1" . $e_value;
            $e_value_list{$e_value} = 1;
        }
    }
    elsif ($input =~ /Sequences producing significant alignments/) { 
        $take_e_value = 1;
    }
    elsif (($take_e_value) and ($input =~ /^>\S+/)) { 
        $take_e_value = 0;
    }
}

foreach $e_value (sort { $a <=> $b } keys %e_value_list) { print "$e_value\n"; }

