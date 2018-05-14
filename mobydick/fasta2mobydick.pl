#!/usr/bin/env perl

# fasta2mobydick.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/21/2008.
# Purpose: convert a FASTA DNA sequence file to a mobydick sequence file.
# Caveat: not at all well-designed to handle huge sequences; use tied files for that.

use strict;
use warnings;

my $seqname    = q{};
my %fasta_seqs = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > (\S+) /xms ) { 
        $seqname = $1;

        # Mobydick file format doesn't allow internal colons:
        $seqname =~ s/:/_/g;

        # Defense against 2+ identical names, e.g., from ": to _" reformat:
        if ( exists $fasta_seqs{$seqname} ) {
            my $i = 1;
            $seqname .= ".$i";
            while ( exists $fasta_seqs{$seqname} ) { 
                $seqname =~ s/\.\d+\z//;
                $i++;
                $seqname .= ".$i";
            }
        }
    }
    elsif ( $input =~ / [a-zA-Z] /xms ) { 
        $input =~ tr/[A-Z]/[a-z]/;
        if ($seqname) { 
            $fasta_seqs{$seqname} .= $input;
        }
    }
}

foreach my $sequence (sort keys %fasta_seqs ) { 
    my $length = length( $fasta_seqs{$sequence} ); 
    print $sequence, 
          ':', 
          $length, 
          ':-', 
          $length, 
          ':', 
          $fasta_seqs{$sequence}, 
          ":\n", 
          ;
}

sub failsafe_name {
    my $filename = $_[0];
    if (-e $filename) {
        my $suffix = 0;
        while (-e $filename) {
            $suffix++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix";
        }
    }
    return $filename;
}

