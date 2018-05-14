#!/usr/bin/env perl

# list_three_TreeFam_group_types_12aug2013.pl -- Erich Schwarz <ems394@cornell.edu>, 8/12/2013.
# Purpose: given hacked_TreeFams/*.newick.txt files extracted from TreeFam v. 9.0, extract their accession numbers, and print out to four groups...
#     human + D. melanogaster + c. elegans ("human_worm.TFaccs.txt");
#     human + non-deuterostome + at least one Caenorhabditis, but no C. elegans ("human_invert_near.miss.worm.TFaccs.txt");
#     human + non-deuterostome + no Caenorhabditis at all ("human_invert_no.worm.TFaccs.txt");
#     human + C. elegans + no D. melanogaster ("human_worm_no.fly.TFaccs.txt").
# Strictly speaking, I only need to distinguish human-worm from human-invert-non-worm, but I am curious to see how many TreeFams are near-misses or no-flies-but-worms.

use strict;
use warnings;

my %c_elegans_match = (
    6239 => 1,   # caenorhabditis_elegans
);

my %human_match = ( 
    9606 => 1,   # homo_sapiens
);

my %melanogaster_match = ( 
    7227 => 1,     # drosophila_melanogaster  
);

my %non_elegans_but_caeno = (
    135651 => 1,   # caenorhabditis_brenneri
    473542 => 1,   # caenorhabditis_briggsae
    281687 => 1,   # caenorhabditis_japonica
    31234 => 1,    # caenorhabditis_remanei
);

# This list was extracted from:  http://www.treefam.org/static/download/treefam_species_tree9.phy, then hacked into NCBI IDs from TreeFam TF*.orthoxml.xml files.
# Note that human, melanogaster, and *all* caenorhabditis have been omitted from this hash.

my %non_worm_non_deuterostome = ( 
    7029 => 1,     # acyrthosiphon_pisum
    7159 => 1,     # aedes_aegypti
    400682 => 1,   # amphimedon_queenslandica
    43151 => 1,    # anopheles_darlingi
    7165 => 1,     # anopheles_gambiae
    7460 => 1,     # apis_mellifera
    3702 => 1,     # arabidopsis_thaliana
    12957 => 1,    # atta_cephalotes
    7091 => 1,     # bombyx_mori
    6326 => 1,     # bursaphelenchus_xylophilus
    283909 => 1,   # capitella_teleta
    7176 => 1,     # culex_quinquefasciatus
    13037 => 1,    # danaus_plexippus
    6669 => 1,     # daphnia_pulex
    7217 => 1,     # drosophila_ananassae
    7220 => 1,     # drosophila_erecta
    7222 => 1,     # drosophila_grimshawi
    7230 => 1,     # drosophila_mojavensis
    7234 => 1,     # drosophila_persimilis
    46245 => 1,    # drosophila_pseudoobscura
    7238 => 1,     # drosophila_sechellia
    7240 => 1,     # drosophila_simulans
    7244 => 1,     # drosophila_virilis
    7260 => 1,     # drosophila_willistoni
    7245 => 1,     # drosophila_yakuba
    34740 => 1,    # heliconius_melpomene
    6412 => 1,     # helobdella_robusta
    37862 => 1,    # heterorhabditis_bacteriophora
    6945 => 1,     # ixodes_scapularis
    225164 => 1,   # lottia_gigantea
    6305 => 1,     # meloidogyne_hapla
    81824 => 1,    # monosiga_brevicollis
    7425 => 1,     # nasonia_vitripennis
    45351 => 1,    # nematostella_vectensis
    121224 => 1,   # pediculus_humanus
    54126 => 1,    # pristionchus_pacificus
    218847 => 1,   # proterospongia_sp
    559292 => 1,   # saccharomyces_cerevisiae
    6183 => 1,     # schistosoma_mansoni
    284812 => 1,   # schizosaccharomyces_pombe
    34506 => 1,    # strongyloides_ratti
    7070 => 1,     # tribolium_castaneum
    6334 => 1,     # trichinella_spiralis
    10228 => 1,    # trichoplax_adhaerens
);

my @input_files = @ARGV;

my $hum_worm       = 'human_worm.cryptic.TFaccs.txt';
$hum_worm          = safename($hum_worm);

my $hum_nearworm   = 'human_invert_near.miss.worm.cryptic.TFaccs.txt';
$hum_nearworm      = safename($hum_nearworm);

my $hum_invert     = 'human_invert_no.worm.cryptic.TFaccs.txt';
$hum_invert        = safename($hum_invert);

my $hum_worm_nofly = 'human_worm_no.fly.cryptic.TFaccs.txt';
$hum_worm_nofly    = safename($hum_worm_nofly);

open my $HUM_WORM,       '>', $hum_worm       or die "Can't open human-worm TF list output $hum_worm: $!";
open my $HUM_NEARWORM,   '>', $hum_nearworm   or die "Can't open human-nearworm TF list output $hum_nearworm: $!";
open my $HUM_INVERT,     '>', $hum_invert     or die "Can't open human-invertebrate TF list output $hum_invert: $!";
open my $HUM_WORM_NOFLY, '>', $hum_worm_nofly or die "Can't open human-worm_no-fly TF list output $hum_worm_nofly: $!";

foreach my $infile (@input_files) {
    if ( $infile =~ /(TF\d+)/xms ) { 
        my $tf_id   = $1;
        my %seen    = ();
        
        open my $INFILE, '<', $infile or die "Can't open input file $infile: $!";
        while (my $input = <$INFILE>) { 
            chomp $input;
            while ( $input =~ / G = [^=]+ T = ( \d+ ) \D /xmsg ) {
                my $species = $1;
                if ( $human_match{$species} ) {
                    $seen{'homo_sapiens'} = 1;
                }
                elsif ( $c_elegans_match{$species} ) {
                    $seen{'elegans'}      = 1;
                    $seen{'invertebrate'} = 1;
                }
                elsif ( $non_elegans_but_caeno{$species} ) {
                    $seen{'caenorhabditis'} = 1;
                    $seen{'invertebrate'}   = 1;
                }
                elsif ( $melanogaster_match{$species} ) {
                    $seen{'melanogaster'} = 1;
                    $seen{'invertebrate'} = 1;
                }
                elsif ( $non_worm_non_deuterostome{$species} ) { 
                    $seen{'invertebrate'} = 1;
                }
            }
        }
        close $INFILE or die "Can't close filehandle to input file $infile: $!";

        # Only print out the TF acc. if both human and some invertebrate are both in the TF group.
        if ( $seen{'homo_sapiens'} and $seen{'invertebrate'} ) { 
            # With C. elegans, two choices:
            if ( $seen{'elegans'} ) { 
                if ( $seen{'melanogaster'} ) {
                    print $HUM_WORM "$tf_id\n";
                }
                else { 
                    print $HUM_WORM_NOFLY "$tf_id\n";
                }
            }
            # Or w/ Caeno. non-elegans:
            elsif ( $seen{'caenorhabditis'} ) { 
                print $HUM_NEARWORM "$tf_id\n";
            }
            # Or plain default:
            else { 
                print $HUM_INVERT "$tf_id\n";
            }
        }
    }
    else { 
        die "Can't extract TreeFam accession from name of input file $infile\n";
    }
}

close $HUM_WORM       or die "Can't close filehandle to human-worm TF list output $hum_worm: $!";
close $HUM_NEARWORM   or die "Can't close filehandle to human-nearworm TF list output $hum_nearworm: $!";
close $HUM_INVERT     or die "Can't close filehandle to human-invertebrate TF list output $hum_invert: $!";
close $HUM_WORM_NOFLY or die "Can't close filehandle to human-worm_no-fly TF list output $hum_worm_nofly: $!";

sub safename {
    my $_filename = $_[0];
    my $_orig_filename = $_filename;
    if (-e $_orig_filename) {
        my $_suffix1 = 1;
        $_filename = $_filename . ".$_suffix1";
        while (-e $_filename) {
            $_suffix1++;
            $_filename =~ s/\.\d+\z//xms;
            $_filename = $_filename . ".$_suffix1";
        }
    }
    return $_filename;
}

