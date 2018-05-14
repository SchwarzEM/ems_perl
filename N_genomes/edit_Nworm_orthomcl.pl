#!/usr/bin/perl

# edit_Nworm_orthomcl.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/2/2009.
# Purpose: turn N-worm protein-based orthomcl.out files to gene-based ones; presupposes WBGene-aware proteome FASTAs.

use strict;
use warnings;

my $input = q{};
my %cds2gene = ();

unless ($#ARGV == 5) { 
    die 'Format: edit_Nworm_orthomcl.pl',
        ' [wormpep, or headers]',          # $ARGV[0] => @proteomes; 
        ' [brigpep, or headers]',          # $ARGV[1]  "
        ' [remapep, or headers]',          # $ARGV[2]  "
        ' [brepep, or headers]',           # $ARGV[3]  "
        ' [jappep, or headers]',           # $ARGV[4]  "
        ' [orthomcl.out]',                 # $ARGV[5] => $ortho_out
        "\n",
        ; 
}

my $ortho_out = pop @ARGV;
my @proteomes = @ARGV;
my @species   = qw( elegans briggsae remanei brenneri japonica );

foreach my $i (0..4) { 
    open my $PROTEOME, '<', $proteomes[$i] 
        or die "Cannot open proteome $proteomes[$i]: $!";
    while (my $input = <$PROTEOME>) { 
        chomp $input;
        if ($input =~ /\A > (\S+) \s+.*\s+ (WBGene\d+) \s+ /xms) { 
            map_cds2gene( $1, $2, $species[$i] );
        }
    }
    close $PROTEOME 
        or die "Can't close filehandle to proteome $proteomes[$i]: $!";
}

# Sample input:
# 
# ORTHOMCL3711(5 genes,4 taxa):    CBP24868(briggsae) Contig48.032(pb2801) F44B9.4a(elegans) F44B9.4b(elegans) cr01.sctg49.wum.86.1(remanei)

open my $ORTHO_OUT, '<', $ortho_out 
    or die "Cannot open N-worm-species OrthoMCL output $ortho_out: $!";
while ($input = <$ORTHO_OUT>) { 
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

        # Zero these out for every single input line:
        my %species_seen = ();
        my %genes_seen = ();

        # Count taxa, and list genes instead of proteins:
        foreach my $o_prot (@orthoprots) { 
            # Fail loudly if either value's missing:
            if ( $o_prot !~ / .+ \( \w+ \) /xms ) {
                die "Can't map $o_prot to species tag!\n";
            }
            if (! $cds2gene{$o_prot} ) {
                die "Can't map $o_prot to gene!\n";
            }
            # If everything's OK, then do the mappings:
            if ( $o_prot =~ / .+ \( (\w+) \) /xms ) { 
                my $species_tag = $1;
                $species_seen{$species_tag} = 1; 
            }
            $genes_seen{$cds2gene{$o_prot}} = 1; 
        }

        # Get summary and final gene list:
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
close $ORTHO_OUT 
    or die "Cannot close filehandle to N-worm-species OrthoMCL output $ortho_out: $!";

sub map_cds2gene { 
        my $_cds_label  = $_[0];
        my $_gene_label = $_[1];   # If equiv., can be set == to $_[0].
        my $_species    = $_[2];

        $_cds_label  .= "($_species)";
        $_gene_label .= "($_species)";

        # This subroutine populates a global data structure, %cds2gene:
        $cds2gene{$_cds_label} = $_gene_label;

        # Pro forma:
        return;
}

