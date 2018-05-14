#!/usr/bin/env perl

# edit_Nworm_orthomcl_cdhit.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/6/2009.
# Purpose: w/ CD-HIT, turn prot.-ish orthomcl.out to nonred.-gene-ish one -- assumes WBGene-aware FASTA headers.

use strict;
use warnings;
use Carp;
use File::Basename;
use Getopt::Long;

my %proteomes      = ();
my %cdhits         = ();
my $orthomcl_input = q{};

# Standard order of recording/reading species data:
my @spp = qw( elegans briggsae remanei brenneri japonica );

my %synonym2cds     = ();
my %cds2gene        = ();

GetOptions (  'ele:s' => \$proteomes{'elegans'},
              'bri:s' => \$proteomes{'briggsae'},
              'rem:s' => \$proteomes{'remanei'},
              'bre:s' => \$proteomes{'brenneri'},
              'jap:s' => \$proteomes{'japonica'},

              'crem:s' => \$cdhits{'remanei'},
              'cbre:s' => \$cdhits{'brenneri'},
              'cjap:s' => \$cdhits{'japonica'},

              'omcl:s' => \$orthomcl_input,         
           );

### Parse input files: ###

my $proteome_count = keys %proteomes;
my $cdhit_count    = keys %cdhits;

unless (     ( $proteome_count >= 2 ) 
         and ( $cdhit_count >= 1   ) 
         and ( $orthomcl_input     ) ) { 
    die_loudly();
}

### Parse CD-HIT data. ###

# Sample input lines:
# 
# >Cluster 88
# 0       2237aa, >Contig1436.japonica.jigsaw.2... *
# >Cluster 89
# 0       424aa, >Contig2391.japonica.jigsaw.3... at 2:414:1319:1731/96%
# 1       2230aa, >Contig95.japonica.jigsaw.17... *

foreach my $cdhit_species_id ( sort keys %cdhits ) { 
    my $key_cds = q{};
    my @other_cdses = ();
    open my $CDHIT, '<', $cdhits{$cdhit_species_id} 
        or croak "Cannot open CD-HIT clustering file",
                 " $cdhits{$cdhit_species_id}: $!",
                 ;

    while (my $input = <$CDHIT>) {
        chomp $input;

        # Commit just-read key and alias CDS names to memory;
        #     get ready to read more.
        if ($input =~ /\A >Cluster \s \d+ \s* /xms) { 

            # Check against failure of CDS mapping:
            if ( (! $key_cds) and (@other_cdses) ) { 
                carp "The alias CDSes \"@other_cdses\" in",
                     " $cdhits{$cdhit_species_id} lack a key",
                     " CDS, at line: \"$input\".\n",
                     ;
            }

            # Map CDS aliases to key CDS name:
            if ($key_cds) { 
                # If any aliases exist, that is...
                if (@other_cdses) {  
                    foreach my $o_c (@other_cdses) {  
                        $synonym2cds{$o_c} = $key_cds;
                    }
                }
                # But key CDS is always its own 'alias'.
                $synonym2cds{$key_cds} = $key_cds;
            }

            # Zero out these values, before reading the next cluster:
            $key_cds = q{};
            @other_cdses = ();
        }

        # Read the key CDS name of a new cluster.
        elsif ( $input =~ / \A 
                            \d+ 
                            .+ 
                            > (\S.+) 
                            [.]{3} \s+ 
                            [*]
                          /xms ) { 

            # Pull in entire line from after '>' up to '... *'.
            $key_cds = $1;

            # Trim off any surplus with greedy "s/\s.\z//;".
            $key_cds =~ s/\s.*\z//;

            # Important! ensure the OrthoMCL suffix:
            if ( $key_cds !~ / \( $cdhit_species_id \) \z /xms ) { 
                $key_cds .= "($cdhit_species_id)";
            }

            # Requires that $key_cds be read only once from all files:
            if ( $synonym2cds{$key_cds} ) { 
                croak "Redundant CDS name in CD-HIT input(s).\n";
            }
        }

        # Read an alias CDS name of a new cluster.
        elsif ( $input =~ / \A 
                            \d+ 
                            .+ 
                            > (\S+) 
                            [.]{3} \s+ 
                            at \s \d 
                          /xms ) { 

            # Again, pull in post-'>' to "... at \d".
            my $o_c = $1;

            # And then trim off any detritus:
            $o_c =~ s/\s.*\z//;

            # Again, ensure OrthoMCL suffixes:
            if ( $o_c !~ / \( $cdhit_species_id \) \z /xms ) { 
                $o_c .= "($cdhit_species_id)";
            }

            # Halt! if any CDS name seen in any previous clusters:
            if ( $synonym2cds{$o_c} ) { 
                croak "Redundant CDS name in CD-HIT input(s).\n";
            }

            # Temporarily array the nonredundant CDS names:
            push @other_cdses, $o_c;
        }

        # Loudly fail if the above parsing fails even once.
        else {
            croak "In $cdhits{$cdhit_species_id}, failed to parse: $input\n";
        }

    }

    # Clear stored data after reading last line of <$CDHIT>:
    if ($key_cds) { 
        if (@other_cdses) {
            foreach my $o_c (@other_cdses) {  
                $synonym2cds{$o_c} = $key_cds;
            }
        }
        $synonym2cds{$key_cds} = $key_cds;
    }
    close $CDHIT 
        or croak "Can't close filehandle to $cdhits{$cdhit_species_id}: $!";
}

### Parse proteomes, using CD-HIT redundancies. ###

foreach my $species (@spp) { 
    if ($proteomes{$species}) { 
        open my $PROTEOME, '<', $proteomes{$species} 
            or croak "Cannot open proteome $proteomes{$species}: $!";
        while (my $input = <$PROTEOME>) {
            chomp $input;
            if ($input =~ /\A > (\S+) \s+.*\s+ (WBGene\d+) \s+ /xms) {
                my $cds    = $1;
                my $wbgene = $2;

                # Ensure OrthoMCL suffix on $cds name:
                if ( $cds !~ / \( $species \) \z  /xms ) { 
                    $cds .= "($species)";
                }

                # And ensure OrthoMCL suffix on $wbgene name:
                if ( $wbgene !~ / \( $species \) \z  /xms ) {
                    $wbgene .= "($species)";
                }

                if ( ( exists $synonym2cds{$cds}         )
                        and ( $synonym2cds{$cds} ne $cds ) ) { 
                    $cds = $synonym2cds{$cds};
                }

                # Make sure that *every* CDS is its own synonym!
                if (     (! $cdhits{$species}   ) 
                     and (! $synonym2cds{$cds}  ) 
                     and ( $proteomes{$species} ) ) {
                    $synonym2cds{$cds} = $cds; 
                }
                map_cds2gene( $cds, $wbgene, $species, \%cds2gene );
            }
            elsif ( $input =~ /\A > /xms ) { 
                croak "Can't parse FASTA header $input!\n";
            }
        }
        close $PROTEOME
            or croak "Can't close filehandle to proteome",
                     " $proteomes{$species}: $!",
                     ;
    }
}

### Parse a single orthomcl.out file, and print out its revision: ###

# Sample input:
# ORTHOMCL5940(5 genes,5 taxa):	 CBG28141(briggsae) CBN08237(brenneri) CJA27759(japonica) CRE00189(remanei) T13C5.8(elegans)

open my $ORTHO_INPUT, '<', $orthomcl_input 
    or croak "Cannot open N-species OrthoMCL output $orthomcl_input: $!";
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
        my @orthoprots = split /\s+/, $oprots_line;
        my %species_seen = ();
        my %genes_seen = ();
        foreach my $o_prot (@orthoprots) { 
            if ( $o_prot !~ / .+ \( \w+ \) /xms ) {
                croak "Can't parse OrthoMCL protein $o_prot\n";
            }
            if ( $o_prot =~ / .+ \( (\w+) \) /xms ) { 
                my $species_tag = $1;
                $species_seen{$species_tag} = 1; 
            }
            $genes_seen{ $cds2gene{ $synonym2cds{$o_prot} } } = 1;
        }
        my $ortho_no  = scalar(keys %genes_seen);
        my $taxon_no  = scalar(keys %species_seen);
        my @orthologs = sort keys %genes_seen;
        print "$ortho_grp",
              "($ortho_no genes,$taxon_no taxa)",
              ":\t",
              " @orthologs\n",
              ;
    }
}
close $ORTHO_INPUT 
    or croak "Cannot close filehandle to $orthomcl_input: $!";


### Subroutines: ###

sub die_loudly {
    die 'Format: ./edit_3worm_orthomcl.pl               \ ', "\n",
        '    # Need at least two proteomes:             \ ', "\n",
        '    --ele  [elegans proteome]                  \ ', "\n",
        '    --bri  [briggsae proteome]                 \ ', "\n",
        '    --rem  [remanei proteome]                  \ ', "\n",
        '    --bre  [brenneri proteome]                 \ ', "\n",
        '    --jap  [japonica proteome]                 \ ', "\n",
        '    # Need at least one CD-HIT *.clustr file:  \ ', "\n",
        '    --crem [remanei *.clustr]                  \ ', "\n",
        '    --cbre [brenneri *.clustr]                 \ ', "\n",
        '    --cjap [japonica *.clustr]                 \ ', "\n",
        '    # Need a single OrthoMCL output:           \ ', "\n",
        '    --ocml [orthomcl.out]                      \ ', "\n",
        "\n",
        ;
}

sub map_cds2gene { 
        my $cds_label        = $_[0];
        my $gene_label       = $_[1];   # If equiv., can be set == to $_[0].
        my $species          = $_[2];
        my $cds2gene_hashref = $_[3];

        # This suffix-labelling should have been done beforehand.
        # So, defend against redundancy.  But, do ensure '($species)' suffix:
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

