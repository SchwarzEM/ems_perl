#!/usr/bin/env perl

# list_three_TreeFam_group_types_12aug2013.pl -- Erich Schwarz <ems394@cornell.edu>, 8/12/2013.
# Purpose: given *.aln.emf files from TreeFam v. 9.0, extract their accession numbers, and print out to three groups...
#     human + D. melanogaster + c. elegans ("human_worm.TFaccs.txt");
#     human + non-deuterostome + at least one Caenorhabditis, but no C. elegans ("human_invert_near.miss.worm.TFaccs.txt");
#     human + non-deuterostome + no Caenorhabditis at all ("human_invert_no.worm.TFaccs.txt");
#     human + C. elegans + no D. melanogaster ("human_worm_no.fly.TFaccs.txt").
# Strictly speaking, I only need to distinguish human-worm from human-invert-non-worm, but I am curious to see how many TreeFams are near-misses or no-flies-but-worms.

use strict;
use warnings;

# This list was extracted from:  http://www.treefam.org/static/download/treefam_species_tree9.phy
# Note that caenorhabditis_elegans has been omitted!
my %non_worm_non_deuterostome = ( 
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
);

my @input_files = @ARGV;

my $hum_worm       = 'human_worm.TFaccs.txt';
$hum_worm          = safename($hum_worm);

my $hum_nearworm   = 'human_invert_near.miss.worm.TFaccs.txt';
$hum_nearworm      = safename($hum_nearworm);

my $hum_invert     = 'human_invert_no.worm.TFaccs.txt';
$hum_invert        = safename($hum_invert);

my $hum_worm_nofly = 'human_worm_no.fly.TFaccs.txt';
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
            if ( $input =~ /\A SEQ \s+ (\S+) \s+/xms ) {
                my $species = $1;
                if ( $species eq 'homo_sapiens' ) {
                    $seen{'homo_sapiens'} = 1;
                }
                elsif ( $species eq 'caenorhabditis_elegans' ) { 
                    $seen{'elegans'}      = 1;
                    $seen{'invertebrate'} = 1;
                }
                elsif ( $species =~ /\Acaenorhabditis_/xms ) { 
                    $seen{'caenorhabditis'} = 1;
                    $seen{'invertebrate'}   = 1;
                }
                elsif ( $species eq 'drosophila_melanogaster' ) { 
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

