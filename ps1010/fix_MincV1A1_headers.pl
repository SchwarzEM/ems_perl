#!/usr/bin/env perl

# fix_MincV1A1_headers.pl -- Erich Schwarz <emsch@its.caltech.edu>, 1/27/2010.
# Purpose: correct the headers of MincV1A1 to have simple clean protein names, and unambiguous 'gene' names.

# Sample input:
# >gnl|MincDB|prot:Minc05158 length:799 contig:MiV1ctg127 region:1761-9463 strand:-

use strict;
use warnings;

my $header    = q{};
my $protname  = q{};
my $contig    = q{};
my $nt_start  = q{};
my $nt_stop   = q{};
my $nt_coords = q{};
my $strand    = q{};
my $gene      = q{};

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > ( gnl\|MincDB\|prot: 
                       (\w+) 
                       \s+ length:\d+ \s+ 
                       contig: (\w+) \s+ region: (\d+) \- (\d+) \s+ strand:([+-]) )
                       \s* \z /xms ) { 
        $header    = $1;
        $protname  = $2;
        $contig    = $3;
        $nt_start  = $4;
        $nt_stop   = $5;
        $strand    = $6;
        if ( $strand eq '+' ) { 
            $nt_coords = $nt_start . q{.} . $nt_stop;
        }
        if ( $strand eq '-' ) {
            $nt_coords = $nt_stop . q{.} . $nt_start;
        }
        $gene      = $contig . '_' . $nt_coords;
        $input     = '>'. "$protname  gene:$gene  $header";
    }
    elsif ( $input =~ / \A > /xms ) { 
        die "Unparseable header:  $input\n";
    }
    print "$input\n";
}

