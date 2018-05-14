#!/usr/bin/env perl

# filter_monster_WB_gff.pl -- Erich Schwarz <emsch@its.caltech.edu>, 1/5/2011.
# Purpose: given WormBase's huge GFFs, get out just the information required to map elements with respect to genes.

use strict;
use warnings;

my $wanted_text = q{};

while (my $input = <>) { 
    chomp $input;
    # Definitely need "Coding_transcript \t protein_coding_primary_transcript", but apprently need lots of other stuff too.
    if ( $input =~ / \A 
                     CHROMOSOME_ 
                     ( \S+ 
                       \t 
                       (?: curated\S* | Coding_transcript | ncRNA ) 
                       \t 
                       (?: exon | coding_exon | intron | CDS | \S*transcript\S* ) 
                       \t 
                       .* 
                     ) \z  /xms ) { 
        $wanted_text = $1;
        print "$wanted_text\n";
    }
    elsif ( $input =~ /\A CHROMOSOME_ 
                          ( \S+ \t gene \t gene \t .* ) 
                       \z /xms ) {
        $wanted_text = $1;
        print "$wanted_text\n";
    }
    # Catch these, too:
    elsif ( $input =~ /\A CHROMOSOME_ 
                          ( \S+ \t [^\t]* \t Pseudogene \t .* \t Pseudogene [ ] \" [^\']+ \" .*) 
                       \z /xms ) { 
        $wanted_text = $1;
        print "$wanted_text\n";
    }
}

