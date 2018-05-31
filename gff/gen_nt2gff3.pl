#!/usr/bin/env perl

# gen_nt2gff3.pl -- Erich Schwarz <ems394@cornell.edu>, 5/31/2018.
# Purpose: convert a DNA sequence file into a GFF3 (useful for single-nucleotide liftovers).

use strict;
use warnings;
use autodie;

use File::Basename;

my $infile   = $ARGV[0];
my $basename = basename($infile);

my $seq = q{};
my $i   = 0;

open my $INFILE, '<', $infile;

while (my $input = <$INFILE>) {
    chomp $input;
    if ( $input =~ /\A > (\S+) /xms ) { 
        $seq = $1;
        $i   = 0;
    }
    elsif ( $input =~ /\S/xms ) {
        $input =~ s/\s//g;
        my @residues = split //, $input;
        foreach my $nt (@residues) {
            $i++;
            print "$seq\t",
                  "$basename\t",
                  "nt\t",
                  "$i\t",
                  "$i\t",
                  ".\t.\t.\t",
                  "seq=$nt",
                  "\n",
                  ;
        }
    }
}

