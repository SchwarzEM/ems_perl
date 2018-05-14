#!/usr/bin/env perl

# extract_anomalous_WashU_gene_coords.pl -- Erich Schwarz <emsch@its.caltech.edu>, 7/27/2011.

# Purpose: map names of anomalous sp. 7/9 genes to coordinates on WS226 scaffolds, given GFF3s.
# Detailed Purpose: given lists like anomalous_Csp7_genes_27jul2011.txt, anomalous_Csp9_genes_27jul2011.txt, 
#     and possibly_anomalous_Csp7_genes_27jul2011.txt, along with the WS226 GFF3s for C. sp. 7 and 9, 
#     get contigs and coordinates for each named anomalous gene.

use strict;
use warnings;

my $gene     = q{};
my $contig   = q{};
my $start_nt = q{};
my $end_nt   = q{};

my $badlist  = $ARGV[0];
my $gff3     = $ARGV[1];
my %badgenes = ();

# List of 'bad' genes:
# g10202

open my $BADLIST, '<', $badlist or die "Can't open list of anomalous or possibly anomalous C. sp. 7/9 WashU genes: $!";
while (my $input = <$BADLIST>) { 
    chomp $input;
    if ( $input =~ /\A (g\d+) \z/xms ) { 
         $gene = $1;
         $badgenes{$gene} = 1;
    }
    else { 
        die "From $badlist, can't parse gene name: $input\n";
    }
}
close $BADLIST or die "Can't close filehandle to list of anomalous or possibly anomalous C. sp. 7/9 WashU genes: $!";

# GFF3 with *lots* of lines, including these ones that I want:
# Contig0 WormBase        gene    17399   22687   .       +       .       ID=gene:g1 ...

open my $GFF3, '<', $gff3 or die "Can't open GFF3 file $gff3: $!";
while (my $input = <$GFF3>) { 
    chomp $input;
    if ( $input =~ /\A ( (?:Contig|Scaffold) \d+) \t WormBase \t gene \t (\d+) \t (\d+) \t \. \t [+|-] \t \. \t ID=gene:(g\d+) /xms ) { 
        $contig   = $1;
        $start_nt = $2;
        $end_nt   = $3;
        $gene     = $4; 
        if ( exists $badgenes{$gene} ) { 
            print "Anomalous gene $gene maps to scaffold $contig, nt $start_nt-$end_nt\n";
            delete $badgenes{$gene};
        }       
    }
}
close $GFF3 or die "Can't close filehandle to GFF3 file $gff3: $!";

foreach my $undetected_gene (sort keys %badgenes) { 
    warn "WARNING -- could not map: $undetected_gene\n";
}

