#!/usr/bin/env perl

# edit_Nworm_orthomcl_cdhit.pl -- Erich Schwarz <emsch@its.caltech.edu>, 7/18/2008.
# Purpose: w/ CD-HIT, turn prot.-ish orthomcl.out to nonred.-gene-ish one.

use strict;
use warnings;
use File::Basename;

my %cdhit_inputs    = ();
my %synonym2cds     = ();
my %cds2gene        = ();
my %proteome_inputs = ();
my $orthomcl_input  = q{};

my %spp_w_cds2gene = ( elegans => 1, );


### Parse @ARGV: ###

unless ($#ARGV >= 2) { 
    die_loudly();
}

foreach my $input_file (@ARGV) { 
    my $filename = basename($input_file);

    # Require '.clstr' at end of name of CD-HIT data.
    if ( $filename =~ / \.clstr \z /xms ) {
        $cdhit_inputs{$input_file} = 1;
    }

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
    if (     ( $filename !~ / \.clstr \z /xms       )
         and ( $filename !~ /orthomcl\.out/xms      )
         and ( $filename !~ /\A ([a-zA-Z]+) \. /xms ) ) { 
        print "Can't parse filename $filename!\n";
        die_loudly();
    }
}


### Halt if any of these three sanity checks are failed: ###

my $cdhit_filecount = scalar (keys %cdhit_inputs);
if ($cdhit_filecount == 0) { 
    die "No CD-HIT *.clstr input files.\n";
}

if (! $orthomcl_input) { 
    die "No OrthoMCL file to change.\n";
}

my @continue = grep { $spp_w_cds2gene{$_}  } 
               map  { $proteome_inputs{$_} }
               keys %proteome_inputs;
if (! @continue) { 
    my @failed_taxa = sort keys %proteome_inputs;
    warn "Only had these taxa: @failed_taxa\n";
    die "No significant changes to make in OrthoMCL.\n";
}


### Parse CD-HIT data. ###

# Sample input lines:
# 
# >Cluster 88
# 0       2237aa, >Contig1436.japonica.jigsaw.2... *
# >Cluster 89
# 0       424aa, >Contig2391.japonica.jigsaw.3... at 2:414:1319:1731/96%
# 1       2230aa, >Contig95.japonica.jigsaw.17... *

foreach my $cdhit_file ( sort keys %cdhit_inputs ) { 
    warn "Reading CD-HIT file $cdhit_file.\n";
    my $key_cds = q{};
    my @other_cdses = ();
    open my $CDHIT, '<', $cdhit_file 
        or die "Cannot open CD-HIT clustering file $cdhit_file: $!";

    while (my $input = <$CDHIT>) {
        chomp $input;
        if ($input =~ /\A >Cluster \s \d+ \s* /xms) { 
            if (! $key_cds) { 
                warn "No key CDS in $cdhit_file before: \"$input\".\n";
            }
            if ($key_cds) {
                if (@other_cdses) {
                    foreach my $o_c (@other_cdses) { 
                        $synonym2cds{$o_c} = $key_cds;
                    }
                }
                # Important!  Or later mapping fails:
                $synonym2cds{$key_cds} = $key_cds;  
            }
            $key_cds = q{};
            @other_cdses = ();
        }
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
            if ( $synonym2cds{$key_cds} ) { 
                die "Redundant CDS name in CD-HIT input(s).\n";
            }
        }
        elsif ( $input =~ / \A 
                            \d+ 
                            .+ 
                            > (\S+.+) 
                            [.]{3} \s+ 
                            at \s \d 
                          /xms ) { 
            # Again, pull in post-'>' to "... at \d".
            my $o_c = $1;
            # And then trim off any detritus:
            $o_c =~ s/\s.*\z//;
            if ( $synonym2cds{$o_c} ) {
                die "Redundant CDS name in CD-HIT input(s).\n";
            }
            push @other_cdses, $o_c;
        }
        elsif ( $input =~ / \A
                            \d+
                            .+
                            > (\S+.+)
                            [.]{3} \s+
                          /xms ) {
            die "In $cdhit_file, failed to parse: $input\n";
        }
        else {
            die "In $cdhit_file, failed to parse: $input\n";
        }

    }

    # Clear stored data after last line of <$CDHIT>
    foreach my $o_c (@other_cdses) { 
        $synonym2cds{$o_c} = "$key_cds";
    }
    $cds2gene{$key_cds} = "$key_cds";
    close $CDHIT or die "Can't close filehandle to $cdhit_file: $!";
}


### Parse proteomes, using CD-HIT redundancies. ###

foreach my $proteome_file ( sort keys %proteome_inputs ) { 
    open my $PROTEOME, '<', $proteome_file
        or die "Cannot open proteome file $proteome_file: $!";
    warn "$proteome_file is assumed to be the",
         " $proteome_inputs{$proteome_file} proteome.\n",
         ;

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
        # But, treat CD-HIT condensations as CDS-to-gene mappings.
        if ( $species ne 'elegans' ) { 
            if ( $input =~ /\A > (\S+) \s* /xms ) { 
                my $cds;
                $cds = $1;
                my $gene = $cds;
                if ( ( exists $synonym2cds{$cds}         ) 
                        and ( $synonym2cds{$cds} ne $cds ) ) {
                    $gene = $synonym2cds{$cds};
                }
                map_cds2gene( $cds, $gene, $species, \%cds2gene );
            }
        }
    }
    close $PROTEOME or die "Can't close filehandle for $proteome_file: $!";
}


### Parse a single orthomcl.out file, and print out its revision: ###

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
            if ( $o_prot !~ / .+ \( \w+ \) /xms ) {
                die "Can't parse OrthoMCL protein $o_prot\n";
            }
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


### Subroutines: ###

sub die_loudly {
    print 'Format: ./edit_3worm_orthomcl.pl           \ ', "\n",
          '    [/optional_directory/species\.\S*]{1+} \ ', "\n",
          '    [/optional_directory/*orthomcl.out*]   \ ', "\n",
          '    [/optional_directory/*.clstr]{1+}     \ ', "\n",
          "\n",
          ;
    die 'Must have *orthomcl.out* and *.clstr in names, \ ', "\n",
        '    but arguments can be in any order.',             "\n",
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

