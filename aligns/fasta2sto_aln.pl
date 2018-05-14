#!/usr/bin/env perl

# fasta2sto_aln.pl -- from Bio::AlignIO perldoc, by Erich Schwarz <emsch@its.caltech.edu>, 2/27/2009.
# Purpose: convert FASTA to Stockholm alignments, for (e.g.) HMMER.

use strict;
use warnings;
use Bio::AlignIO;
use File::Basename;

foreach my $input_align (@ARGV) { 
    my $output_align = basename($input_align);
    $output_align =~ s/\.fa\z//;
    $output_align .= '.sto';

    if (-r $input_align) { 
        if (-e $output_align ) { 
            die "Won't overwrite $output_align!\n";
        }
        my $in  = Bio::AlignIO->new( -file   => "$input_align",
                                     -format => 'fasta',        );
        my $out = Bio::AlignIO->new( -file => ">$output_align", 
                                     -format => 'stockholm',    );
        while ( my $aln = $in->next_aln() ) { 
            $out->write_aln($aln);
        }
    }
}

