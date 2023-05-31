#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = '##gff-version 3';

while (my $input = <>) {
    chomp $input;
    # Sample input:
    # chrI      2       82      81      6       AGGCTT
    if ( $input =~ /\A (\S+) \t (\d+) \t (\d+) \t (\d+) \t (\d+) \t ([ACGT]+) \z/xms ) {
        # Definitions from: https://www.ensembl.org/info/website/upload/gff3.html 
        # Fields must be tab-separated. 
        # Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.'
        # More GFF3 format documenation at: http://gmod.org/wiki/GFF3
        # The SOFA sequence ontology is browsable at: http://www.sequenceontology.org/browser/obob.cgi
        # 'tandem_repeat' is defined at: http://www.sequenceontology.org/browser/current_release/term/SO:0000705

        my $seqid       = $1;               # name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix.
        my $source      = 'mTR';            # name of the program that generated this feature, or the data source (database or project name)
        my $type        = 'tandem_repeat';  # type of feature. Must be a term or accession from the SOFA sequence ontology
        my $start       = $2;               # Start position of the feature, with sequence numbering starting at 1.
        my $end         = $3;               # End position of the feature, with sequence numbering starting at 1.
        my $score       = q{.};             # A floating point value.
        my $strand      = q{+};             # defined as + (forward) or - (reverse).
        my $phase       = q{.};             # One of '0', '1' or '2' [if defined at all; otherwise '.'].
        my $tr_length   = $4;
        my $unit_length = $5;
        my $unit_seq    = $6;

        # attributes - A semicolon-separated list of tag-value pairs, providing additional information about each feature. 
        # Some of these tags are predefined, e.g. ID, Name, Alias, Parent - see the GFF documentation for more details.
        my $attributes  = "tr_length=$tr_length; unit_length=$unit_length; unit_seq=$unit_seq;";

        print "$header\n" if $header;
        $header = q{};

        print "$seqid\t$source\t$type\t$start\t$end\t$score\t$strand\t$phase\t$attributes\n";
    }
    else {
        die "Cannot parse input: $input\n";
    }
}
