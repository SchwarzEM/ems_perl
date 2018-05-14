#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $name    = q{};
my $prefix  = q{};
my $type    = q{};
my $species = q{};
my $source  = q{};

my %rep2sp = (
    Agam => 'Anopheles gambiae',
    Amel => 'Apis mellifera',
    Asuu => 'Ascaris suum',
    Bmal => 'Brugia malayi',
    Bmor => 'Bombyx mori',
    Btau => 'Bos taurus',
    Bxyl => 'Bursaphelenchus xylophilus',
    Cang => 'Caenorhabditis angaria',
    Cbri => 'Caenorhabditis briggsae',
    Cele => 'Caenorhabditis elegans', 
    Cfam => 'Canis familiaris',
    Cgig => 'Crassostrea gigas',
    Cint => 'Ciona intestinalis',
    Ctel => 'Capitella teleta',
    Dimm => 'Dirofilaria immitis',
    Drer => 'Danio rerio',
    Ggal => 'Gallus gallus',
    Ggor => 'Gorilla gorilla',
    Hbac => 'Heterorhabditis bacteriophora',
    Hcon => 'Haemonchus contortus',
    Hsap => 'Homo sapiens',
    Lgig => 'Lottia gigantea',
    Lloa => 'Loa loa',
    Mdel => 'Monodelphis domestica',
    Mhap => 'Meloidogyne hapla',
    Mmul => 'Macaca mulatta',
    Mmus => 'Mus musculus',
    Nam  => 'Necator americanus',
    Oana => 'Ornithorhynchus anatinus',
    Pabe => 'Pongo abelii',
    Pexs => 'Pristionchus exspectatus',
    Phum => 'Pediculus humanus',
    Pmar => 'Petromyzon marinus',
    Ppac => 'Pristionchus pacificus',
    Pred => 'Panagrellus redivivus',
    Ptro => 'Pan troglodytes',
    Sman => 'Schistosoma mansoni',
    Spur => 'Strongylocentrotus purpuratus',
    Tspi => 'Trichinella spiralis',
    Tsui => 'Trichuris suis',
);

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        if ( $input =~ /\A > (\S+ \.RB_ \S+) \t ([^\t]+) \t ([^\t]+) \z/xms ) {
            $name    = $1;
            $type    = $2;
            $species = $3;
            $source  = 'RepBase';
            $type    =~ s/\A\s+//;
            $type    =~ s/\s+\z//;
            $species =~ s/\A\s+//;
            $species =~ s/\s+\z//;
        }
        # Acey_RScout_0001
        elsif ( $input =~ /\A > (Acey_RScout_\d+) \b /xms ) { 
            $name    = $1;
            $type    = 'ab initio prediction';
            $species = 'Ancylostoma ceylanicum';
            $source  = 'This study';
        }
        elsif ( $input =~ /\A > ( (\S+) _R=\d+) \b /xms ) {
            $name    = $1;
            $prefix  = $2;
            $type    = 'ab initio prediction';  
            $species = $rep2sp{$prefix};
            $source  = 'This study';
        }
        else { 
            die "Cannot parse sequence line: $input\n";
        }
        print "$name\t$species\t$type\t$source\n";
    }
}

