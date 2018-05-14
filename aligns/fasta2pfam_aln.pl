#!/usr/bin/env perl

# fasta2pfam_aln.pl -- copied from Bio::AlignIO perldoc by Erich Schwarz <emsch@its.caltech.edu>, 2/27/2009.
# Purpose: convert FASTA to PFAM (Stockholm) alignments, for (e.g.) HMMER.

use strict;
use warnings;
use Bio::AlignIO;
use File::Basename;

foreach my $input_fasta (@ARGV) { 
    my $output_pfam = basename($input_fasta);
    $output_pfam =~ s/\.fa\z//;
    $output_pfam .= '.pfam';

    if (-r $input_fasta) { 
        if (-e $output_pfam ) { 
            die "Won't overwrite $output_pfam!\n";
        }
        my $in  = Bio::AlignIO->new( -file   => "$input_fasta",
                                     -format => 'fasta',        );
        my $out = Bio::AlignIO->new( -file => ">$output_pfam", 
                                     -format => 'pfam',         );
        while ( my $aln = $in->next_aln() ) { 
            $out->write_aln($aln);
        }
    }
}

