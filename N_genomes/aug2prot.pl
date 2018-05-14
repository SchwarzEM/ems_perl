#!/usr/bin/env perl

# aug2prot.pl -- Erich Schwarz <emsch@its.caltech.edu>, 4/15/2009.
# Purpose: get FASTA of protein translations from *.aug report.

use strict;
use warnings;

my $transcript= q{};
my $readprot  = 0;
my %name2seq = ();

while (my $input = <>) { 
    chomp $input;
    my $residues;
    # Sample input: [...] transcript_id "NODE_34_length_4412_cov_14.754.g1.t2"; [...]
    # Note: this gets read again and again, line after line;
    #    the idea is to capture the *last* line before the protein seq. starts.
    if ($input =~ / transcript_id \s+ \"([^\"]+)\" /xms) { 
        $transcript = $1;
    }
    # Sample input: # protein sequence = [GQR ...  KQK]  <-- all seq. on one line
    if ( (! $readprot) and ( $input =~ / \A
                                         [#]{1} \s+
                                         protein \s+
                                         sequence \s+
                             = \s+
                             \[
                             ([A-Z]+)
                             \] \s* \z /xms ) ) {
        $residues = $1;
        $name2seq{$transcript} .= $residues;
        $readprot = 0;
    }
    # Sample input: # protein sequence = [MSK...LHF  <-- seq. only *starts* on one line
    if ( (! $readprot) and ( $input =~ / \A 
                                         [#]{1} \s+ 
                                         protein \s+ 
                                         sequence \s+ 
                             = \s+ 
                             \[ 
                             ([A-Z]+) 
                             \s* \z /xms ) ) { 
        $residues = $1;
        $name2seq{$transcript} .= $residues;
        $readprot = 1;
    }
    # Sample input: # NHL...DRK]
    if ( ($readprot) and ( $input =~ / \A [#]{1} \s+ 
                                     ([A-Z]+)
                                     \] \s* \z /xms ) ) { 
        $residues = $1;
        $name2seq{$transcript} .= $residues;
        $readprot = 0;
    }
    # Sample input: # PRR...EQI
    if ( ($readprot) and ( $input =~ / \A [#]{1} \s+ 
                                     ([A-Z]+) 
                                     \s* \z /xms ) ) { 
        $residues = $1;
        $name2seq{$transcript} .= $residues;
    }
}

foreach my $protein (sort keys %name2seq) { 
    print ">$protein\n";
    my @output_lines 
        = unpack("a60" x (length($name2seq{$protein})/60 + 1), $name2seq{$protein});
    foreach my $output_line (@output_lines) { 
        if ($output_line =~ /\S/) { 
            print "$output_line\n";
        }
    }
}

