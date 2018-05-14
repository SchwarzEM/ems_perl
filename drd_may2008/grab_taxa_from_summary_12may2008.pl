#!/usr/bin/env perl

# grab_taxa_from_summary_12may2008.pl -- Erich Schwarz <emsch@its.caltech.edu>
# Purpose: a deliberately hand-written script aimed at selecting taxa of "biological interest" from a homolog summary.

use strict;
use warnings;

# N.B. many Bacterial taxa were omitted for the sake of brevity, which should go into a more detailed study.
# Bacterial choices were heavily driven by authorial bias -- what E.M.S. remembered hearing something about.

my @terms = ( 
              'Methanosarcina mazei Go1',               # Archaea; Euryarchaeota; Methanomicrobia
              'Methanosphaera stadtmanae DSM 3091',     # Archaea; Euryarchaeota; Methanobacteria

              'Acidobacteria bacterium Ellin345',            # Bacteria; Acidobacteria; Acidobacteriales; Acidobacteriaceae
              'Streptomyces coelicolor A3(2)',               # Bacteria; Actinobacteria; Actinobacteridae; Actinomycetales
              'Bacteroides thetaiotaomicron VPI-5482',       # Bacteria; Bacteroidetes; Bacteroidetes (class); Bacteroidales
              'Chlorobium chlorochromatii CaD3',             # Bacteria; Chlorobi; Chlorobia; Chlorobiales
              'Chloroflexus aurantiacus J-10-fl',            # Bacteria; Chloroflexi; Chloroflexales
              'Synechococcus sp. CC9311',                    # Bacteria; Cyanobacteria; Chroococcales; Synechococcus
              'Deinococcus radiodurans R1',                  # Bacteria; Deinococcus-Thermus; Deinococci; Deinococcales
              'Bacillus subtilis subsp. subtilis str. 168',  # Bacteria; Firmicutes; Bacillales; Bacillaceae
              'Mycoplasma penetrans HF-2',                   # Bacteria; Firmicutes; Mollicutes; Mycoplasmataceae
    'Fusobacterium nucleatum subsp. polymorphum ATCC 10953', # Bacteria; Fusobacteria; Fusobacteriales; Fusobacteriaceae
              'Rhodopirellula baltica SH 1',                 # Bacteria; Planctomycetes; Planctomycetacia; Planctomycetales
              'Rhizobium leguminosarum bv. viciae 3841',     # Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales
              'Rickettsia typhi str. Wilmington',            # Bacteria; Proteobacteria; Alphaproteobacteria; Rickettsiales
              'Neisseria gonorrhoeae FA 1090',               # Bacteria; Proteobacteria; Betaproteobacteria; Neisseriales
              'Myxococcus xanthus DK 1622',                  # Bacteria; Proteobacteria; Deltaproteobacteria; Myxococcales
              'Escherichia coli str. K12',                   # Bacteria; Proteobacteria; Gammaproteobacteria
              'Pseudomonas aeruginosa PAO1',                 # Bacteria; Proteobacteria; Gammaproteobacteria; Pseudomonadales; Pseudomonadaceae; Pseudomonas; Pseudomonas aeruginosa PAO1
              'Arcobacter butzleri RM4018',                  # Bacteria; Proteobacteria; Epsilonproteobacteria
              'Magnetococcus sp. MC-1',                      # Bacteria; Proteobacteria; Magnetococcus
            'Leptospira interrogans serovar Lai str. 56601', # Bacteria; Spirochaetes; Spirochaetales; 

              'Cryptosporidium hominis TU502',      # Eukaryota; Alveolata; Apicomplexa
              'Monosiga brevicollis MX1',           # Eukaryota; Choanoflagellida
              'Aspergillus oryzae RIB40',           # Eukaryota; Fungi; Dikarya; Ascomycota; Pezizomycotina; Eurotiomycetes
              'Sclerotinia sclerotiorum 1980',      # Eukaryota; Fungi; Dikarya; Ascomycota; Pezizomycotina; Leotiomycetes
              'Neurospora crassa OR74A',            # Eukaryota; Fungi; Dikarya; Ascomycota; Pezizomycotina; Sordariomycetes
              'Ustilago maydis 521',                # Eukaryota; Fungi; Dikarya; Basidiomycota; Ustilaginomycotina
              'Mus musculus',                       # Eukaryota; Metazoa; Chordata; Craniata; Vertebrata
              'Caenorhabditis elegans',             # Eukaryota; Metazoa; Nematoda
              'Drosophila melanogaster',            # Eukaryota; Metazoa; Arthropoda
              'Nematostella vectensis',             # Eukaryota; Metazoa; Cnidaria; Anthozoa
              'Dictyostelium discoideum AX4',       # Eukaryota; Mycetozoa; Dictyosteliida
              'Ostreococcus lucimarinus CCE9901',   # Eukaryota; Viridiplantae; Chlorophyta

              'Enterobacteria phage Sf6',           # Viruses; dsDNA viruses, no RNA stage; Caudovirales; Podoviridae; P22-like viruses
              'Pseudomonas phage D3',               # Viruses; dsDNA viruses, no RNA stage; Caudovirales; Siphoviridae
            );
# Perl Cookbook, 2cd. ed., pp. 203-205:
my @patts = map { qr/$_/ } @terms;

while (my $input = <>) { 
    chomp $input;
    SCAN:
    foreach my $pattern (@patts) { 
        if ($input =~ /$pattern/xms) { 
            # Any match? print once, and move on:
            print "$input\n";
            last SCAN;
        }
    }
}

