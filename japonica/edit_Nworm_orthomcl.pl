#!/usr/bin/env perl

# edit_Nworm_orthomcl.pl -- Erich Schwarz <emsch@its.caltech.edu>, 7/18/2008.
# Purpose: turn protein-based orthomcl.out files to gene-based ones.

use strict;
use warnings;
use File::Basename;

my %cds2gene        = ();
my %proteome_inputs = ();
my $orthomcl_input  = q{};

my %significant_spp = ( elegans => 1 );

unless ($#ARGV >= 2) { 
    die_loudly();
}

foreach my $input_file (@ARGV) { 
    my $filename = basename($input_file);

    # Require 'orthomcl.out' within name of OMCL output.
    if ( $filename =~ /orthomcl\.out/xms ) { 
        $orthomcl_input = $input_file;
    }

    # Typical proteome filename: 'elegans.fa'.
    if (     ( $filename !~ /orthomcl\.out/xms      ) 
         and ( $filename =~ /\A ([a-zA-Z]+) \. /xms ) ) { 
        my $species = $1;
        $proteome_inputs{$input_file} = $species;
    }

    # Failsafe:
    if (     ( $filename !~ /orthomcl\.out/xms      )
         and ( $filename !~ /\A ([a-zA-Z]+) \. /xms ) ) { 
        print "Can't parse filename $filename!\n";
        die_loudly();
    }
}

# Halt if these sanity checks are failed:
if (! $orthomcl_input) { 
    die "No OrthoMCL file to change!\n";
}
my @continue = ();
@continue = grep { $significant_spp{$_} } 
            map { $proteome_inputs{$_} }
            keys %proteome_inputs;
if (! @continue) { 
    my @failed_taxa = sort keys %proteome_inputs;
    warn "Only had these taxa: @failed_taxa\n";
    die "No significant changes to make in OrthoMCL.\n";
}

foreach my $proteome_file ( sort keys %proteome_inputs ) { 
    open my $PROTEOME, '<', $proteome_file
        or die "Cannot open proteome file $proteome_file: $!";
    warn "$proteome_file is assumed to be the $proteome_inputs{$proteome_file} proteome.\n";

    while (my $input = <$PROTEOME>) { 
        chomp $input;
        my $species = $proteome_inputs{$proteome_file};
        if ( $species eq 'elegans' ) { 
            if ( $input =~ /\A > (\S+) \s+.*\s+ (WBGene\d+) \s+ /xms ) { 
                my ($cds, $wbgene);
                $cds    = $1;
                $wbgene = $2;
                map_cds2gene( $cds, $wbgene, $species, \%cds2gene );
            }
        }
        # Default: gene:protein 1:1 ratio; 'protein name' eq 'gene name'.
        if ( $species ne 'elegans' ) { 
            if ( $input =~ /\A > (\S+) \s* /xms ) { 
                my $cds;
                $cds = $1;
                map_cds2gene( $cds, $cds, $species, \%cds2gene );
            }
        }
    }
    close $PROTEOME or die "Can't close filehandle for $proteome_file: $!";
}


# Sample input:
# 
# ORTHOMCL3711(5 genes,4 taxa):    CBP24868(briggsae) Contig48.032(pb2801) F44B9.4a(elegans) F44B9.4b(elegans) cr01.sctg49.wum.86.1(remanei)

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
        my @orthoprots = split /\s+/, $oprots_line;
        my %species_seen = ();
        my %genes_seen = ();
        foreach my $o_prot (@orthoprots) { 
            if ( $o_prot =~ / .+ \( (\w+) \) /xms ) { 
                my $species_tag = $1;
                $species_seen{$species_tag} = 1; 
            }
            $genes_seen{$cds2gene{$o_prot}} = 1;
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
close $ORTHO_INPUT or die "Cannot close filehandle to $orthomcl_input: $!";

sub die_loudly {
    die 'Format: ./edit_3worm_orthomcl.pl',
        ' [/optional_directory/species\.\S*]{1+ separate files}',
        ' [*orthomcl.out* (must have that in name, but argument order not crucial)]',
        "\n",
        ;
}

sub map_cds2gene { 
        my $cds_label        = $_[0];
        my $gene_label       = $_[1];   # If equiv., can be set == to $_[0].
        my $species          = $_[2];
        my $cds2gene_hashref = $_[3];

        $cds_label  .= "($species)";
        $gene_label .= "($species)";

        # This subroutine populates *any* hash passed in as a reference.
        # Originally, it populated a global data structure, %cds2gene.
        $cds2gene_hashref->{$cds_label} = $gene_label;

        # Pro forma:
        return;
}

