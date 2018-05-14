#!/usr/bin/env perl

# make_gff3_header_pragma.pl -- Erich Schwarz <ems@emstech.org>, 2/10/2014.
# Purpose: given the FASTA sequence file that a GFF3 is meant to annotate, generate a GFF3 header text with pragmas for all sequences in that text.  Based on http://www.sequenceontology.org/gff3.shtml, date 26 February 2013, version 1.21.

use strict;
use warnings;

my $seq = q{};

my $data_ref;

my $header = '##gff-version 3';

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > (\S+) /xms ) { 
        $seq = $1;
    }
    else { 
        $input =~ s/\s//g;
        $data_ref->{'seq'}->{$seq}->{'residues'} .= $input;
    }
}

my @sequences = sort keys %{ $data_ref->{'seq'} };

foreach my $seq (@sequences) { 
     print "$header\n" if $header;
     $header = q{};
     my $seq_length = 0;
     if (! exists $data_ref->{'seq'}->{$seq}->{'residues'} ) {
         die "Can't identify (or count) residues of sequence $seq\n";
     }
     if ( exists $data_ref->{'seq'}->{$seq}->{'residues'} ) {
         my $residues = $data_ref->{'seq'}->{$seq}->{'residues'};
         $seq_length = length($residues);
     }
     print '##sequence-region ', "$seq 1 $seq_length\n";
}

