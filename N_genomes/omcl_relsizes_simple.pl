#!/usr/bin/env perl

# omcl_relsizes_simple.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/17/2009.
# Purpose: for one or more OrthoMCLs, simple table w/ ratios of each sp. longest protein to elegans's longest protein.
# Note: rough-and-ready way to see if one sibling species' gene predictions are more fragmentary than the others.

use strict;
use warnings;
use Getopt::Long;

my %proteomes      = ();
my %prot2len       = ();
my $orthomcl_input = q{};
my $show_abs_nos;

# Standard order of recording/reading species data:
my @spp = qw( elegans briggsae remanei brenneri japonica );

GetOptions (  'ele:s' => \$proteomes{'elegans'},
              'bri:s' => \$proteomes{'briggsae'},
              'rem:s' => \$proteomes{'remanei'},
              'bre:s' => \$proteomes{'brenneri'},
              'jap:s' => \$proteomes{'japonica'},

              'omcl:s' => \$orthomcl_input,         
              'abs'    => \$show_abs_nos,
           );

### Parse input files: ###

my $proteome_count = keys %proteomes;

unless (     ( $proteome_count >= 2  ) 
         and ( $proteomes{'elegans'} )
         and ( $orthomcl_input       ) ) { 
    die_loudly();
}

### Read and store protein lengths for each proteome. ### 

foreach my $species (@spp) { 
    if ($proteomes{$species}) { 
        open my $PROTEOME, '<', $proteomes{$species} 
            or die "Cannot open proteome $proteomes{$species}: $!";
        my $cds = q{};
        while (my $input = <$PROTEOME>) {
            chomp $input;
            if ($input =~ /\A > (\S+) \s* /xms) {
                $cds    = $1;
                # Ensure OrthoMCL suffix on $cds name:
                if ( $cds !~ / \( $species \) \z  /xms ) { 
                    $cds .= "($species)";
                }
                # And then start the residue count.
                $prot2len{$cds} = 0;
            }
            elsif ( $input =~ /\S/xms ) {
                # Absolutely minimal filtering:
                $input =~ s/\s//g;
                # Increment residue count:
                $prot2len{$cds} += length($input);
            }
        }
        close $PROTEOME
            or die "Can't close filehandle to proteome",
                   " $proteomes{$species}: $!",
                   ;
    }
}

### Given orthomcl.out file lines, get max-res counts and ratios vs. elegans: ###

# Sample input:
# ORTHOMCL5940(5 genes,5 taxa):	 CBG28141(briggsae) CBN08237(brenneri) CJA27759(japonica) CRE00189(remanei) T13C5.8(elegans)

open my $ORTHO_INPUT, '<', $orthomcl_input 
    or die "Cannot open N-species OrthoMCL output $orthomcl_input: $!";

while (my $input = <$ORTHO_INPUT>) { 
    chomp $input;
    if ($input =~ / \A
                    (ORTHOMCL\d+)                 # $1 -> $ortho_grp
                    \(\d+\sgenes,\d+\staxa\)      # just (punctuation)
                    : \s+ 
                    (.+)                          # $2 -> $oprots_line
                    \s* 
                    \z 
                  /xms) { 
        my $ortho_grp = $1;
        my $oprots_line = $2;
        my %species2max = ();
        foreach my $sp1 (@spp) { 
            $species2max{$sp1} = 0;
        }
        my @orthoprots = split /\s+/, $oprots_line;
        foreach my $o_prot (@orthoprots) { 
            if ( $o_prot !~ / .+ \( \w+ \) /xms ) {
                die "Can't parse OrthoMCL protein $o_prot\n";
            }
            if ( $o_prot =~ / .+ \( (\w+) \) /xms ) { 
                my $species_tag = $1;
                if (! $prot2len{$o_prot} ) { 
                    die "Failed to record length of protein $o_prot!\n";
                }
                my $prot_size = $prot2len{$o_prot};
                if ( $prot_size > $species2max{$species_tag} ) { 
                    $species2max{$species_tag} = $prot_size;
                }
            }
        }
        my @length_data = ();
        push @length_data, $ortho_grp;
        if ( $species2max{'elegans'} >= 1 ) { 
            foreach my $sp2 (@spp) { 
                my $ratio = ($species2max{$sp2} / $species2max{'elegans'} );
                $ratio = sprintf "%.2f", $ratio;
                my $datum = "$ratio";
                if ($show_abs_nos) { 
                    $datum .= " [$species2max{$sp2}]";
                }
                push @length_data, $datum;
            }
            my $output = join "\t", @length_data;
            print "$output\n";
        }
    }
}
close $ORTHO_INPUT 
    or die "Cannot close filehandle to $orthomcl_input: $!";


### Subroutines: ###

sub die_loudly {
    die 'Format: ./edit_3worm_orthomcl.pl                 \ ', "\n",
        '    # Need elegans and at least other proteome:  \ ', "\n",
        '    --ele  [elegans proteome, REQUIRED]          \ ', "\n",
        '    --bri  [briggsae proteome]                   \ ', "\n",
        '    --rem  [remanei proteome]                    \ ', "\n",
        '    --bre  [brenneri proteome]                   \ ', "\n",
        '    --jap  [japonica proteome]                   \ ', "\n",
        '    # Need a single OrthoMCL output:             \ ', "\n",
        '    --ocml [orthomcl.out]                        \ ', "\n",
        '    # Optional flag to insert max. aa no.:       \ ', "\n",
        '    --abs  [show absolute aa. nos.]              \ ', "\n",
        ;
}

sub map_cds2gene { 
        my $cds_label        = $_[0];
        my $gene_label       = $_[1];   # If equiv., can be set == to $_[0].
        my $species          = $_[2];
        my $cds2gene_hashref = $_[3];

        # Defend against redundancy, but always append '(species)':
        if ( $cds_label !~ / \( $species \) \z/xms ) { 
            $cds_label  .= "($species)";
        } 
        if ( $gene_label !~ / \( $species \) \z/xms ) { 
            $gene_label .= "($species)";
        }

        # This subroutine populates *any* hash passed in as a reference.
        # Originally, it populated a global data structure, %cds2gene.
        $cds2gene_hashref->{$cds_label} = $gene_label;

        # Pro forma:
        return;
}

