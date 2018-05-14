#!/usr/bin/env perl

# rework_ngasp_gffseqs.pl -- Erich Schwarz <emsch@its.caltech.edu>, 7/11/2008.
# Purpose: rework nGASP-assoc. GFFs into Sequence or CDS/gene_id lines that I can BioPerl.
# 
# N.B.: this silently deletes 'CDS' lines that span entire genes 
#     rather than being exons.  That probably matters!

use strict;
use warnings;

$^I = ".bak";

my $contig = q{};
my $bulk_line = q{};

# N.B.: this silently deletes 'CDS' lines that span entire genes rather than being exons.
# That may have been a factor in getting the reformatting to work!

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ / \A \s* \# /xms) { 
        print "$input\n";
    }

    # Pass on Sequence length-lines unchanged.
    if (  ( $input =~ / \A
                        ( [^\t]+ ) 
                        \t \. \t
                        Sequence
                        \t 1 \t \d+ \t .+ 
                        \z 
                      /xms ) 
          and  ( $1 =~ /\S/ ) ) { 
        print "$input\n";
    }

    # For briggsae.jigsaw.gff2:
    # Alter Chr\d+ / CDS / Transcript lines to 
    #    Chr\d+ CDS/gene_id, print; silently ignore most others.
    # 
    # Sample input:
    # chrI nGASP CDS 28321 28380 . -  Transcript ContigchrI.briggsae.jigsaw.3.mRNA3

    if ( $input =~ / \A    
                   ( ( [^\t]+ ) \t [^\t]* \t CDS \t (?: [^\t]* \t){5} ) # $part1 ($seqname)
                   Transcript                                           # to 'gene_id'
                   \s+ 
                   ( \S+ )                                              # $CDS_name
                   /xms ) {

        my ($part1, $seqname, $CDS_name);
        $part1    = $1; 
        $seqname  = $2;
        $CDS_name = $3;  
        if ( $seqname !~ /\S/ ) {
            die "Unusable sequence name in:\n$input\n";
        }
        $input = $part1 . 'gene_id "'. $CDS_name . q{"} ;
        print "$input\n";
    }              

    # Alter exon/CDS lines to CDS/gene_id, print; silently ignore most others.
    if ( $input =~ / \A
                   ( Contig\d+ \t [^\t]* \t )  # to $part1
                   exon                        # to 'CDS'
                   ( \t (?: [^\t]* \t){5} )    # to $part2
                   CDS                         # to 'gene_id'
                   \s+ \" 
                   ( [^\"]+ ) 
                   \" \s*    # to $part3
                   /xms ) {
        my ($part1, $part2, $CDS_name);
        $part1     = $1;
        $part2     = $2;
        $CDS_name = $3;
        $input = $part1 
                 . 'CDS' 
                 . $part2 
                 . "gene_id \"$CDS_name\""
                 ;
        print "$input\n";
    }
}

