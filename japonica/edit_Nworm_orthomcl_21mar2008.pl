#!/usr/bin/perl

# edit_Nworm_orthomcl.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/21/2008.
# Purpose: turn N-worm protein-based orthomcl.out files to gene-based ones.

use strict;
use warnings;

my $input = "";
my %cds2gene = ();

# unless ($#ARGV == 5) { 
#     die 'Format: edit_3worm_orthomcl.pl',
#         ' [wormpep, or headers]',            # $ARGV[0] -> $WORMPEP
#         ' [brigpep2gid, or headers]',        # $ARGV[1] -> $BRIGPEP2GID
#         ' [rempep, or headers]',             # $ARGV[2] -> $REMPEP
#         ' [brenpep, or headers]',            # $ARGV[3] -> $BRENPEP
#         ' [japonpep, or headers];,           # $ARGV[4] -> $JAPONPEP
#         ' [orthomcl.out]',                   # $ARGV[5] -> $ORTHO_OUT
#         "\n",
#         ; 
# }

open my $WORMPEP, "<", "$ARGV[0]" 
    or die "Cannot open wormpep $ARGV[0]: $!";
while ($input = <$WORMPEP>) { 
    chomp $input;
    if ($input =~ /\A > (\S+) \s+.*\s+ (WBGene\d+) \s+ /xms) { 
        map_cds2gene( $1, $2, 'elegans' );
    }
}
close $WORMPEP;

# Sample brigpep2gid header-line input:
#
# >CBP00000    WBGene00023535; CBG00020    TR:A8WM56    protein_id:CAP21560.1

open my $BRIGPEP2GID, '<', "$ARGV[1]" 
    or die "Cannot open brigpep2gid $ARGV[1]: $!";
while ($input = <$BRIGPEP2GID>) { 
    chomp $input;
    if ($input =~ / \A > (CBP\d+) \s+ (WBGene\d+) /xms) { 
        map_cds2gene( $1, $2, 'briggsae' );
    }
}
close $BRIGPEP2GID;

# Sample input:
# 
# >cr01.sctg0.wum.1004.1 Contig0f-snap.94.final
# >cr01.sctg0.wum.1004.2 Contig0f.Fgenesh_Celegans.63.final

open my $REMPEP, "<", "$ARGV[2]"  
    or die "Cannot open preliminary rempep $ARGV[2]: $!";
while ($input = <$REMPEP>) {
    chomp $input;
    if ($input =~ / \A 
                    > 
                    ( (\S+) \.\d+ )   # CDSes have .\d+ suffixes.
                    \s+ 
                    /xms) { 
        map_cds2gene( $1, $2, 'remanei' );
    }
}
close $REMPEP;

open my $BRENPEP, "<", "$ARGV[3]"  
    or die "Cannot open preliminary brenpep $ARGV[3]: $!";
while ($input = <$BRENPEP>) {
    chomp $input;
    if ($input =~ / \A 
                    > 
                   ( Contig\d+\.\d+\w+ )   # \w+ sp-suffix added, 2008.
                   \s* 
                   /xms) { 
        map_cds2gene( $1, $1, 'brenneri' );
    }
}
close $BRENPEP;

open my $JAPONPEP, "<", "$ARGV[4]" 
    or die "Cannot open preliminary japonpep $ARGV[4]: $!";
while ($input = <$JAPONPEP>) {
    chomp $input;
    if ($input =~ / \A 
                    > 
                    ( Contig\d+\.\d+\w+ )   # \w+ sp-suffix added, 2008.
                    \s* 
                  /xms) { 
        map_cds2gene( $1, $1, 'japonica' );
    }
}
close $BRENPEP;

# Sample input:
# 
# ORTHOMCL3711(5 genes,4 taxa):    CBP24868(briggsae) Contig48.032(pb2801) F44B9.4a(elegans) F44B9.4b(elegans) cr01.sctg49.wum.86.1(remanei)

open my $ORTHO_OUT, "<", "$ARGV[5]" 
    or die "Cannot open N-worm-species OrthoMCL output $ARGV[5]: $!";
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
close $ORTHO_OUT;

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

