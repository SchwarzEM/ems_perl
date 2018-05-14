#!/usr/bin/env perl

# hack_NCBI_taxa_from_TreeFam_XML_18aug2013s.pl -- Erich Schwarz <ems394@cornell.edu>, 8/18/2013.
# Purpose: quickly extract NCBI taxonomy IDs from TreeFam XML, for just the species that I happen to need in a format usable for hashing in Perl.

use strict;  
use warnings;

my %sp2ncbi      = ();
my %seen_species = ();

my %acceptable_species = ( 
    arabidopsis_thaliana => 1,
    saccharomyces_cerevisiae => 1, 
    schizosaccharomyces_pombe => 1,
    saccharomyces_cerevisiae_S288c => 1,   
    'schizosaccharomyces_pombe_972h-' => 1,
    danio_rerio => 1,
    mus_musculus => 1,
    homo_sapiens => 1,
    caenorhabditis_elegans => 1,
    drosophila_melanogaster => 1,
);

while (my $input = <>) {
    chomp $input;
    # Capture this: <species NCBITaxId="7234" name="drosophila_persimilis">
    while ($input =~ /species [ ] NCBITaxId = \" (\d+) \" [ ] name = \" ([^\"\s]+) \"/xmsg ) { 
        my $ncbi_taxon = $1;
        my $species    = $2;
        if ( $acceptable_species{$species} ) { 
            if ( $sp2ncbi{$species} and ( $sp2ncbi{$species} != $ncbi_taxon ) ) {
                die "For species $species, two inconsistent NCBI taxon IDs: $sp2ncbi{$species} versus $ncbi_taxon\n";
            }
            $sp2ncbi{$species}      = $ncbi_taxon;
            $seen_species{$species} = 1;
        }
    }
}

my @final_species = sort keys %sp2ncbi;
foreach my $final_sp (@final_species) {
    print "    $sp2ncbi{$final_sp}", ' => 1,   #', " $final_sp\n";
}

my @possible_species = sort keys %acceptable_species;
foreach my $possible_sp (@possible_species) { 
    if (! $seen_species{$possible_sp} ) { 
        print "NOTE -- failed to detect any instances of species: $possible_sp\n";
    }
}

