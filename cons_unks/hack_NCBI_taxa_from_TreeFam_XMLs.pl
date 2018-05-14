#!/usr/bin/env perl

# hack_NCBI_taxa_from_TreeFam_XMLs.pl -- Erich Schwarz <ems394@cornell.edu>, 8/15/2013.
# Purpose: quickly extract NCBI taxonomy IDs from TreeFam XML, for just the species that I happen to need in a format usable for hashing in Perl.

use strict;  
use warnings;

my %sp2ncbi      = ();
my %seen_species = ();

my %acceptable_species = ( 
    homo_sapiens => 1,
    caenorhabditis_elegans => 1,,
    caenorhabditis_briggsae_AF16 => 1,
    saccharomyces_cerevisiae_S288c => 1,
    'schizosaccharomyces_pombe_972h-' => 1,
    lottia_gigantea => 1,
    capitella_teleta => 1,
    helobdella_robusta => 1,
    ixodes_scapularis => 1,
    atta_cephalotes => 1,
    apis_mellifera => 1,
    nasonia_vitripennis => 1,
    drosophila_virilis => 1,
    drosophila_mojavensis => 1,
    drosophila_grimshawi => 1,
    drosophila_willistoni => 1,
    drosophila_pseudoobscura_pseudoobscura => 1,
    drosophila_persimilis => 1,
    drosophila_yakuba => 1,
    drosophila_simulans => 1,
    drosophila_sechellia => 1,
    drosophila_melanogaster => 1,
    drosophila_erecta => 1,
    drosophila_ananassae => 1,
    anopheles_darlingi => 1,
    anopheles_gambiae => 1,
    culex_quinquefasciatus => 1,
    aedes_aegypti => 1,
    heliconius_melpomene => 1,
    danaus_plexippus => 1,
    bombyx_mori => 1,
    tribolium_castaneum => 1,
    pediculus_humanus_corporis => 1,
    acyrthosiphon_pisum => 1,
    daphnia_pulex => 1,
    trichinella_spiralis => 1,
    pristionchus_pacificus => 1,
    bursaphelenchus_xylophilus => 1,
    meloidogyne_hapla => 1,
    strongyloides_ratti => 1,
    heterorhabditis_bacteriophora => 1,
    caenorhabditis_briggsae => 1,
    caenorhabditis_japonica => 1,
    caenorhabditis_brenneri => 1,
    caenorhabditis_remanei => 1,
    schistosoma_mansoni => 1,
    nematostella_vectensis => 1,
    amphimedon_queenslandica => 1,
    trichoplax_adhaerens => 1,
    saccharomyces_cerevisiae => 1,
    schizosaccharomyces_pombe => 1,
    proterospongia => 1,
    monosiga_brevicollis => 1,
    arabidopsis_thaliana => 1,
    drosophila_pseudoobscura => 1,
    pediculus_humanus => 1,
    proterospongia_sp => 1,
);

# N.B.: added on the last three after an earlier version of this script didn't properly detect them.

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

