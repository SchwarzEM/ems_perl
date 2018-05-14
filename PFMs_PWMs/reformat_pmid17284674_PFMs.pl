#!/usr/bin/env perl

# reformat_Zhao_PFMs.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/31/2010.
# Purpose: given a supplemental data file with PFMs, reformat it into a usable proto-.ace PFM file.

use strict;
use warnings;

my $matrix_no  = q{};
my $consensus  = q{};
my $residue    = q{};
my $values     = q{};
my $site_count = 0;
my $print_head = 1;

# [WBPaper00029109; pmid17284674]

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A Matrix: (\d+) /xms ) { 
        $matrix_no  = $1;
        $consensus  = q{};
        $residue    = q{};
        $values     = q{};
        $site_count = 0;
        if ($print_head) { 
            print "\n";
            $print_head = 0;
        }
        print "Position_Matrix : \"Zhao_et.al_2007_$matrix_no\"\n";
        print "Description \"Matrix no. $matrix_no correlated with muscle-specific gene expression, from Zhao et al. (2007).\"  Paper_evidence \"WBPaper00029109\" \\\\ pmid17284674\n";
        print "Type           Frequency\n";
    }
    elsif ( $input =~ /\A ( [A-Z] (?: \s [A-Z])+ ) \s* \z/xms ) { 
        $consensus = $1;
        $consensus =~ s/\s//g;
    }
    elsif ( $input =~ /\A ( [ACGT] ) \s+ \| \s+ ( (\d+) (?: \s+ \d+ )+ ) \s* \z/xms ) { 
        $residue = $1;
        $values  = $2;
        $site_count += $3;
        print "Site_values   $residue  $values\n";
        if ($residue eq 'T') { 
            print "Sites_used    $site_count\n";
            print "Remark        \"Predicted in WBPaper00029109/pmid17284674; consensus is $consensus.\"\n";
            print "\n";
        }
    }
}

