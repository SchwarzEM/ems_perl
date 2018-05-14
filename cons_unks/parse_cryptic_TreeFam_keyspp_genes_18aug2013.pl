#!/usr/bin/env perl

# parse_covert_TreeFam_keyspp_genes_18aug2013.pl -- Erich Schwarz <ems394@cornell.edu>, 8/18/2013.
# Purpose: given COVERT TF*.newick.txt hacked from TreeFam v. 9.0, make a simple table with columns: TreeFam ID \t key_spp \t gene.  To be then fed to other scripts.

use strict;
use warnings;

my %ok_spp = ( 
    3702   => 'arabidopsis_thaliana',
    6239   => 'caenorhabditis_elegans',
    7955   => 'danio_rerio',
    7227   => 'drosophila_melanogaster',
    9606   => 'homo_sapiens',
    10090  => 'mus_musculus',
    559292 => 'saccharomyces_cerevisiae',
    284812 => 'schizosaccharomyces_pombe',
);

my @input_files = @ARGV;

foreach my $infile (@input_files) {
    if ( $infile =~ /(TF\d+)/xms ) { 
        my $tf_id   = $1;
        open my $INFILE, '<', $infile or die "Can't open input file $infile: $!";
        while (my $input = <$INFILE>) { 
            chomp $input;
            while ( $input =~ / G = ([^=]+) : T = ( \d+ ) \D /xmsg ) {
                my $gene    = $1;
                my $species = $2;
                if ( $ok_spp{$species} ) {
                    # Don't want just to detect wanted species, but rename it from NCBI tax. no. to human-readable form.
                    print "$tf_id\t$ok_spp{$species}\t$gene\n";
                }
            }
        }
        close $INFILE or die "Can't close filehandle to input file $infile: $!";
    }
    else {
        die "Can't extract TreeFam accession from name of input file $infile\n";
    }
}           

