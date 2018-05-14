#!/usr/bin/env perl

# add_pfm_chrome.pl -- Erich Schwarz <emsch@its.caltech.edu>, 6/19/2008.
# Purpose: a simple script for making absolutely minimal ace-like text a bit more usable.

use strict;
use warnings;

my $description = q{Description  "Insert description here."  Paper_evidence "WBPaper00024505" // pmid15492775 };
my $type        = q{Type         Frequency};
my $remark      = q{Remark       "PFM from Supporting Information Dataset S1 of Gaudet et al. (2004) [WBPaper00024505/pmid15492775]."};

my $first_entry = 'yes';
my $in_a_block = 'no';

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A Position_Matrix /xms ) { 
        $in_a_block = 'yes';
        if ($first_entry eq 'yes') { 
            $first_entry = 'no';
            print "\n";
        }
        print "$input\n";
        print "$description\n";
        print "$type\n";
    }
    elsif ( $input =~ /\A ([acgtACGT]) (\s+\d.*) \z /xms ) { 
        my $residue = $1;
        my $line    = $2;
        $residue =~ tr/[a-z]/[A-Z]/;
        print "Site_values  $residue  $line\n";
    }
    elsif ( $input =~ /\S/xms ) { 
        print "$input\n";
    }
    elsif ( ( $input !~ /\S/xms ) and ( $in_a_block eq 'yes' ) ) { 
        $in_a_block = 'no';
        print "$remark\n";
        print "\n";
    }
}

if ( ( eof(ARGV) ) and ( $in_a_block eq 'yes' ) ) {
    print "$remark\n";
    print "\n";
}

