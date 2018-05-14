#!/usr/bin/env perl

# prot2gene_omcl.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/25/2014.
# Purpose: read protein-based OrthoMCL output, then gene table from prot2gene_table.pl, to make gene-based OrthoMCL output.
# 
# Note: this was a modernization of my code from 2009.  It makes the script fail explicitly at the exact point that a CDS-to-gene connection fails,
#    which previously I had not had the sense to make it do.  It also makes failures to parse happen absolutely automatically -- no clumsy double parsing attempts;
#    either a parse works or there is an else-die.
#    It does change the output in some way (probably sorting) that makes exact drop-in replacement of the older version impossible.
# 
# Note: this explicitly gives up any idea of getting easy CDS-to-gene maps from FASTA headers!
#     Making those maps becomes an explicit, prior responsibility of the user.
#     On the other hand, that also means that a separate script for making tables
#         can be devised for each case, archived, and perhaps generalized at some point.
#     It also means that CD-HIT merges, etc., can be captured as de facto CDS-to-gene
#         tables *before* trying to parse raw (protein-based) OrthoMCL outputs.

use strict;
use warnings;
use autodie;
use Getopt::Long;
use List::MoreUtils qw(uniq);

my $input_omcl    = q{};
my @omcl_groups   = ();
my $cds2gene_file = q{};
my $permit_1g1s;
my $data_ref;
my $help;

GetOptions ( 'input_omcl:s' => \$input_omcl,
             'cds2gene=s'   => \$cds2gene_file,
             'permit'       => \$permit_1g1s,
             'help'         => \$help,         );

if ( ($help) or (! -r $input_omcl ) ) { 
    die "Format: prot2gene_omcl.pl\n",
        "        --input_omcl|-i  [output file from OrthoMCL, named \"all_orthomcl.out\" by default]\n",
        "        --cds2gene|-c    [CDS-to-gene table]\n",
        "        --permit|-p      [override default censorship of \"(1 genes,1 taxa)\" orthology groups]\n",
        ;
}

open my $INPUT_OMCL, '<', $input_omcl;
while (my $input = <$INPUT_OMCL>) { 
    chomp $input;
    if ($input =~ / \A
                    (ORTHOMCL\d+)                 # $1 -> $ortho_grp
                    \( \d+\sgenes,\d+\staxa \)    # just (punctuation)
                    : \s+ 
                    (\S.+\S)                      # $2 -> $oprots_line
                    \s* 
                    \z 
                  /xms) { 
        my $ortho_grp   = $1;
        my $oprots_line = $2;

        # For later use, store exact order of groups in parent file:
        push @omcl_groups, $ortho_grp;

        my @orthoprots = split /\s+/, $oprots_line;

        # Archive the orthology groups, as follows:
        # 
        # To enforce only one species per CDS:
        # $data_ref->{'CDS'}->{$cds_id}->{'species'} = $species_tag;
        # 
        # To allow (grudgingly) one or more orthology groups per CDS:
        # $data_ref->{'CDS'}->{$cds_id}->{'omcl_grp'}->{$ortho_grp} = 1;
        # 
        # To list all CDSes belonging to an orthology group, under 'species' to enable per-species listing:
        # $data_ref->{'omcl_grp'}->{$ortho_grp}->{'species'}->{$species_tag}->{'CDS'}->{$cds_id} = 1;
        # 
        # A lethal flag, deleted when CDS-to-gene data is read in:
        # $data_ref->{'no_gene'}->{$cds_id} = 1;

        foreach my $o_prot (@orthoprots) { 
            # If format's OK, extract, check and store various mappings:
            if ( $o_prot =~ / \A ( [^\s\(\)]+ ) \( ( [^\s\(\)]+ ) \) \z /xms ) { 
                my $cds_id     = $1;
                my $species_tag = $2;

                # Enforce the rule: only one species per protein ID!
                if ( exists $data_ref->{'CDS'}->{$cds_id}->{'species'} ) {
                    my @other_spec_list 
                        = grep { $_ ne $species_tag } 
                        sort keys %{ $data_ref->{'CDS'}->{$cds_id}->{'species'} };
                    if (@other_spec_list) {
                        my $spec_list = join ', ', @other_spec_list;
                        die "CDS $cds_id is listed not only for species",
                            " $species_tag, but for the species: $spec_list\n",
                            ;
                    }
                }

                # But in well-formed input, any protein should be seen exactly once.  So:
                if (! exists $data_ref->{'CDS'}->{$cds_id}->{'species'} ) { 
                    # Map protein ID to species:
                    $data_ref->{'CDS'}->{$cds_id}->{'species'} = $species_tag;

                    # It's generally a bad sign -- but not a show-stopper -- if a protein belongs 
                    #     to more than one orthology group.  So emit a warning:
                    if ( exists $data_ref->{'CDS'}->{$cds_id}->{'omcl_grp'} ) { 
                        my @extra_omcl_groups 
                            = sort keys %{ $data_ref->{'CDS'}->{$cds_id}->{'omcl_grp'} };
                        my $extras_list = join ', ', @extra_omcl_groups;
                        warn "CDS $cds_id, in addition to orthology group $ortho_grp,",
                             " belongs to the orthology group(s) $extras_list!\n",
                             ;
                    }

                    # For protein ID, add the orthology group to protein's list
                    #     (ideally, a 1-member list):
                    $data_ref->{'CDS'}->{$cds_id}->{'omcl_grp'}->{$ortho_grp} = 1;

                    # Finally, add new protein (via its species) to orthology group's list of members:
                    $data_ref->{'omcl_grp'}->{$ortho_grp}->{'species'}->{$species_tag}->{'CDS'}->{$cds_id} = 1;

                    # And mark a lethal flag, enforcing a later CDS-to-gene map:
                    $data_ref->{'no_gene'}->{$cds_id} = 1;
                }
            }
            # Enforce correct format:
            else {
                die "Can't parse protein $o_prot in orthology group $ortho_grp of input OrthoMCL file $input_omcl!\n";
            }
        }
    }
    # Enforce correct format:
    else { 
        die "From OrthoMCL file $input_omcl, can't parse input line: $input\n";
    }
}
close $INPUT_OMCL;

# Archive the CDS-to-gene mappings, as follows:
# 
# delete $data_ref->{'no_gene'}->{$cds_id};
# $data_ref->{'CDS'}->{$cds_id}->{'gene'} = $gene;

open my $CDS2GENE, '<', $cds2gene_file;
while (my $input = <$CDS2GENE>) { 
    chomp $input;
    if ($input =~ /\A (\S+) \t (\S+) \z /xms) { 
        my $cds_id  = $1;
        my $gene = $2;
        delete $data_ref->{'no_gene'}->{$cds_id};

        # Note: we're mapping CDSes to genes, *not* proteins to genes.
        #     A single protein *can* be produced by 2+ (very similar) genes!
        #     Mercifully, there really is just one CDS per gene...

        if ($data_ref->{'CDS'}->{$cds_id}->{'gene'}) { 
            my $alt_gene = $data_ref->{'CDS'}->{$cds_id}->{'gene'};
            if ($alt_gene ne $gene) { 
                die "CDS $cds_id assigned to two genes: $alt_gene and $gene!\n";
            }
        }
        $data_ref->{'CDS'}->{$cds_id}->{'gene'} = $gene;
    }
    # Enforce correct format:
    else {
        die "From CDS-to-gene table $cds2gene_file, malformatted input line: $input\n";
    }
}
close $CDS2GENE;

foreach my $ortho_grp (@omcl_groups) { 
    my $output_1    = "$ortho_grp";
    my $output_2    = q{};
    my $output_3    = q{};
    my @o_spp       = sort keys %{ $data_ref->{'omcl_grp'}->{$ortho_grp}->{'species'} };
    my $taxon_count = @o_spp;
    my $gene_count  = 0;

    foreach my $o_species (@o_spp) { 
        my @o_genes = ();
        # First, get the CDSes belonging to the orthology group:
        my @o_cdses = sort keys %{ $data_ref->{'omcl_grp'}->{$ortho_grp}->{'species'}->{$o_species}->{'CDS'} };
        if (! @o_cdses) {
            die "Not able to identify CDSes for OrthoMCL group $ortho_grp in species $o_species\n";
        }

        # It isn't enough to just map from CDSes to genes; we need to make the resulting
        #     list nonredudant, because >= 2 CDSes will often map to 2 copies of a single 
        #     gene name.  These will be stupidly overcounted if not made a unique list.
        foreach my $o_cds (@o_cdses) { 
            my $o_gene = q{};

            if ( exists $data_ref->{'CDS'}->{$o_cds}->{'gene'} ) {
                $o_gene = $data_ref->{'CDS'}->{$o_cds}->{'gene'};
            }
            else { 
                die "For CDS $o_cds from OrthoMCL group $ortho_grp and species $o_species, cannot identify gene\n";
            }
            if (! $o_species) { 
                die "Not able to identify species $o_species for OrthoMCL group $ortho_grp and CDS $o_cds\n";
            }
            $o_gene = $o_gene . q{(} . $o_species . q{)};
            push @o_genes, $o_gene;
        }

        @o_genes    = sort @o_genes;
        @o_genes    = uniq @o_genes;        
        $gene_count += @o_genes;

        # Trying to do direct '.=' to $output_3 fails. And, always use 'q{ }'!
        my $o_gene_line = join q{ }, @o_genes;
        $output_3 .= q{ } . $o_gene_line;
    }
    $output_2 = "($gene_count genes,$taxon_count taxa):\t ";
    if ( ( $output_2 ne "(1 genes,1 taxa):\t " ) or $permit_1g1s ) { 
        print $output_1, $output_2, $output_3, "\n";
    }
}

