#!/usr/bin/perl -w

# rnai_tmaker2go.pl
# Erich Schwarz <emsch@its.caltech.edu>, 3/10/04
# Purpose: convert WS120+ Tablemaker output to gene_association.wb subset.

use strict;

my $input_line = "";
chomp (my @input_lines = <>);

foreach $input_line (@input_lines) 
{
    # de-quote all parts of the input line
    $input_line =~ s/\"//g;

    my @line_part = split /\t/, $input_line;

    # de-suffix the a-z from cosmid.number names for genes
    $line_part[3] =~ s/[a-z]$//;

    # column 1 with obligatory 'WB' database name
    print "WB\t";

    # column 2 with protein name, without leading 'WP:'
    $line_part[4] =~ s/^WP://;
    print "$line_part[4]\t";

    # column 3 with cosmid.N (w/o a-z) or xyz-N name

    if ($line_part[5] =~ /^\s*$/)
    {
        # cosmid.N (w/o a-z) name
        print "$line_part[3]\t";
    }
    else
    {
        # xyz-N name
        print "$line_part[5]\t";
    }

    # column 4 with default ""
    print "\t";

    # column 5 with de-quoted GO:x term
    print "$line_part[1]\t";

    # column 6 with or without reference
    if ($line_part[6]) 
    {
        # with ref.
        print "WB:";
        print "$line_part[6]\t";
    }
    else 
    {
        # or, w/o ref.
        print "\t";
    }

    # column 7 with automatic 'IMP' entry
    print "IMP\t";

    # column 8 with RNAi evidence
    print "WB:";
    print "$line_part[2]\t";

    # column 9 with automatic 'P' for 'biological process'
    print "P\t";
    # note, lazy kludge -- should check GO terms individually for their source ontology

    # column 10 automatically empty
    print "\t";

    # column 11 with gene synonym, if any
    if ($line_part[5] =~ /^\s*$/) 
    {
        print "\t";
    }
    else 
    {
        print "$line_part[3]\t";
    }

    # column 12 with automatic 'gene'
    print "gene\t";

    # column 13 with automatic NCBI taxon code for C. elegans
    print "taxon:6239\t";

    # column 14 with eight-all-digit date
    chomp(my $date = `date +%Y%m%d`);
    print "$date\t";

    # column 15 with highly redundant automatic 'WB', and *end* of line
    print "WB\n"
}
