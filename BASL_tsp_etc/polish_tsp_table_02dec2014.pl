#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %pref2spp = (
    ANOGA => 'Anopheles gambiae',
    ARATH => 'Arabidopsis thaliana',
    BOVIN => 'Bos taurus',
    BRAFL => 'Branchiostoma floridae',
    CAEEL => 'Caenorhabditis elegans',
    CANFA => 'Canis lupus familiaris',
    CHICK => 'Gallus gallus',
    CIOIN => 'Ciona intestinalis',
    CRYNJ => 'Cryptococcus neoformans var. neoformans JEC21',
    DANRE => 'Danio rerio',
    DICDI => 'Dictyostelium discoideum',
    DROME => 'Drosophila melanogaster',
    HUMAN => 'Homo sapiens',
    IXOSC => 'Ixodes scapularis',
    LEIMA => 'Leishmania major',
    MACMU => 'Macaca mulatta',
    MONBE => 'Monosiga brevicollis',
    MONDO => 'Monodelphis domestica',
    MOUSE => 'Mus musculus',
    NEMVE => 'Nematostella vectensis',
    ORNAN => 'Ornithorhynchus anatinus',
    PANTR => 'Pan troglodytes',
    PHYPA => 'Physcomitrella patens',
    RAT   => 'Rattus norvegicus',
    SCHMA => 'Schistosoma mansoni',
    TAKRU => 'Takifugu rubripes',
    TRIVA => 'Trichomonas vaginalis',
    XENTR => 'Xenopus (Silurana) tropicalis',
);

while (my $input = <>) {
    chomp $input;
    # Sample input:
    # TSN11_HUMAN_16-245      sp|A1L157|TSN11_HUMAN/16-245
    if ( $input =~ /\A ( \S+ _ (\S+) _ \d+ [-] \d+) \t (\S+) \/ (\d+ [-] \d+) \z/xms ) {
        my $seqname = $1;
        my $sp_pref = $2;
        my $info    = $3;
        my $range   = $4;
        if (! exists $pref2spp{$sp_pref} ) {
            die "Can't identify species for: $input\n";
        }
        my $species = $pref2spp{$sp_pref};
        print "$seqname\t$species\t$info\t$range\n";
    }
    else {
        die "Can't parse: $input\n";
    }
}

