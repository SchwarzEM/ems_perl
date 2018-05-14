#!/usr/bin/env perl

# slice_GFF_crude.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/17/2009.
# Purpose: very quick, crude script to get GFF2/3 slices.

use strict;
use warnings;

my $arg_no = @ARGV;
if ($arg_no != 4) { 
    die "./slice_GFF_crude.pl [chr] [nt1] [nt2] [GFF] -- if nt1, nt2 both 0, print whole chromosome\n";
}

my %ok_chrs = ( I   => 1,
                II  => 1,
                III => 1,
                IV  => 1,
                V   => 1,
                X   => 1, );
my $chr = $ARGV[0];
if ( (! $chr) or (! $ok_chrs{$chr} ) ) { 
    die "\"$chr\" not OK chromosome\n";
}

my $nt1    = $ARGV[1];
my $nt2    = $ARGV[2];
my $infile = $ARGV[3];

open my $IN, '<', $infile or die "Can't open $infile, wha?\n";
while (my $input = <$IN>) { 
    if ( $input =~ / \A (\S+) \t [^\t]* \t [^\t]* \t (\d+) \t (\d+) /xms ) { 
        my $Nchr = $1;
        my $Nnt1 = $2;
        my $Nnt2 = $3;
        if ( $Nchr eq $chr ) { 
            if ( ( $nt1 == 0 ) and ( $nt2 == 0 ) ) { 
                print $input;
            }
            elsif ( ( $nt1 <= $Nnt1 ) and ( $Nnt2 <= $nt2 ) ) { 
                print $input;
            }
        }
    }
}        

