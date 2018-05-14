#!/usr/bin/env perl

# parse_overt_TreeFam_keyspp_genes_18aug2013.pl -- Erich Schwarz <ems394@cornell.edu>, 8/18/2013.
# Purpose: given *.aln.emf files from TreeFam v. 9.0, make a simple table with columns: TreeFam ID \t key_spp \t gene.  To be then fed to other scripts.

use strict;
use warnings;

my %ok_spp = ( 
    arabidopsis_thaliana => 1,
    saccharomyces_cerevisiae => 1, 
    schizosaccharomyces_pombe => 1,
    danio_rerio => 1,
    mus_musculus => 1,
    homo_sapiens => 1,
    drosophila_melanogaster => 1,
    caenorhabditis_elegans => 1,
);

my @input_files = @ARGV;

foreach my $infile (@input_files) {
    if ( $infile =~ /(TF\d+)/xms ) { 
        my $tf_id   = $1;
        open my $INFILE, '<', $infile or die "Can't open input file $infile: $!";
        while (my $input = <$INFILE>) { 
            chomp $input;
            if ( $input =~ /\A SEQ \s+ (\S+) \s+ (?: \S+ \s+){5} (\S+) /xms ) {
                my $species = $1;
                my $gene    = $2;
                if ( $ok_spp{$species } ) { 
                    print "$tf_id\t$species\t$gene\n";
                }
            }
        }
        close $INFILE or die "Can't close filehandle to input file $infile: $!";
    }
    else { 
        die "Can't extract TreeFam accession from name of input file $infile\n";
    }
}

