#!/usr/bin/env perl

# prot2gene_ofind.pl -- Erich Schwarz <ems394@cornell.edu>, 9/10/2016.
# Purpose: read protein-based OrthoFinder output that *tries* to be OrthoMCL-like (but, in my opinion, fails), then gene table from prot2gene_table.pl, to make gene-based OrthoMCL output.
#
# A separate script for making protein-to-gene tables can be can be devised for each case, if needed; their results can all be pooled before use.
# Also, CD-HIT merges, etc., can be captured as de facto CDS-to-gene tables *before* trying to parse raw (protein-based) OrthoMCL outputs.

use strict;
use warnings;
use autodie;
use Getopt::Long;
use List::MoreUtils qw(uniq);

my $input_ofnd   = q{};
my @ofnd_groups  = ();
my $cds2gene2spp = ();
my $permit_1g1s;
my $data_ref;
my $help;

GetOptions ( 'input_ofnd=s'   => \$input_ofnd,
             'cds2gene2spp=s' => \$cds2gene2spp,
             'permit'         => \$permit_1g1s,
             'help'           => \$help,
);

if ( $help or (! -r $input_ofnd ) or (! -r $cds2gene2spp ) ) {
    die "Format: prot2gene_ofind.pl\n",
        "        --input_ofnd|-i     [output file from OrthoFinder]\n",
        "        --cds2gene2spp|-c   [CDS-to-gene-to-species tab-delimited table; each CDS must be unique]\n",
        "        --permit|-p         [override default censorship of \"(1 genes,1 taxa)\" orthology groups]\n",
        "        --help|-h           [print this message]\n",
        ;
}

open my $CDS2GENE2SPP, '<', $cds2gene2spp;
while (my $input = <$CDS2GENE2SPP>) {
    chomp $input;
    if ($input =~ /\A (\S+) \t (\S+) \t (\S+) \z/xms) {
        my $cds_id  = $1;
        my $gene    = $2;
        my $species = $3;
     
        # Enforce unique CDS names throughout all the species sets being archived.
        # This will also enforce unique mappings of CDS to gene to species.
        #
        # Note: we're mapping CDSes to genes, *not* proteins to genes.
        #     A single protein *can* be produced by 2+ (very similar) genes!
        #     Mercifully, there really is at least one CDS per gene...

        # It turns out that OrthoFinder silently changes ':', '(', and ')' in CDS/protein names to '_';
        #     so this must also be done with CDS/protein names in cds2gene2species.
        $cds_id =~ s/[:]/_/g;
        $cds_id =~ s/[(]/_/g;
        $cds_id =~ s/[)]/_/g;

        if ( exists $data_ref->{'CDS'}->{$cds_id} ) {
            die "Redundant CDS ID: $cds_id\n";
        }

        $data_ref->{'CDS'}->{$cds_id}->{'gene'}    = $gene;
        $data_ref->{'CDS'}->{$cds_id}->{'species'} = $species;

        if ( ( exists $data_ref->{'gene'}->{$gene}->{'species'} ) and ( $species ne $data_ref->{'gene'}->{$gene}->{'species'} ) ) {
            die "Gene $gene mapped to two different species: $species and $data_ref->{'gene'}->{$gene}->{'species'}\n";
        }
        $data_ref->{'gene'}->{$gene}->{'species'} = $species;
    }               
    else {
        die "From CDS-to-gene-to-species tab-delimited table $cds2gene2spp, malformatted input line: $input\n";
    }
}
close $CDS2GENE2SPP;

open my $INPUT_OFND, '<', $input_ofnd;
while (my $input = <$INPUT_OFND>) { 
    chomp $input;
    if ( $input =~ / \A (OG\d+): \s+ (\S.+\S) \s* \z /xms ) { 
        my $ortho_grp   = $1;
        my $oprots_line = $2;

        # For later use, store exact order of groups in parent file:
        push @ofnd_groups, $ortho_grp;

        my @orthoprots = split /\s+/, $oprots_line;

        foreach my $o_prot (@orthoprots) { 
            # Although this should be superfluous for OrthoFinder, to be safe, revise ':', '(', and ')' to '_" anyway.
            $o_prot =~ s/[:]/_/g;
            $o_prot =~ s/[(]/_/g;
            $o_prot =~ s/[)]/_/g;

            if (! exists $data_ref->{'CDS'}->{$o_prot}->{'gene'} ) {
                die "In OrthoFinder group $ortho_grp, cannot map CDS $o_prot to gene.\n";
            }
            my $gene = $data_ref->{'CDS'}->{$o_prot}->{'gene'};

            if (! exists $data_ref->{'CDS'}->{$o_prot}->{'species'} ) {
                die "In OrthoFinder group $ortho_grp, cannot map CDS $o_prot to species.\n";
            }
            my $species = $data_ref->{'CDS'}->{$o_prot}->{'species'};

            # For gene, add the orthology group to protein's list (ideally, a 1-member list):
            $data_ref->{'gene'}->{$gene}->{'ofnd_grp'}->{$ortho_grp} = 1;

            # Finally, add gene (via its species) to orthology group's list of members:
            $data_ref->{'ofnd_grp'}->{$ortho_grp}->{'species'}->{$species}->{'gene'}->{$gene} = 1;
        }
    }
    else { 
        die "From OrthoFinder file $input_ofnd, can't parse input line: $input\n";
    }
}
close $INPUT_OFND;

# Check for genes that belong to two or more OrthoFinder groups.
my @all_genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene (@all_genes) {
    my @orth_grps      = sort keys %{ $data_ref->{'gene'}->{$gene}->{'ofnd_grp'} };
    my $orth_grp_count = @orth_grps;
    if ( $orth_grp_count >= 2 ) {
        my $species = $data_ref->{'gene'}->{$gene}->{'species'};
        warn "MULTIGROUP: Gene $gene from species $species belongs to $orth_grp_count OrthoFinder groups: @orth_grps\n";
    }
}

foreach my $ortho_grp (@ofnd_groups) { 
    my $output_1    = "$ortho_grp";
    my $output_2    = q{};
    my $output_3    = q{};
    my @o_spp       = sort keys %{ $data_ref->{'ofnd_grp'}->{$ortho_grp}->{'species'} };
    my $taxon_count = @o_spp;
    my $gene_count  = 0;

    foreach my $o_species (@o_spp) { 
        my @orig_o_genes = sort keys %{ $data_ref->{'ofnd_grp'}->{$ortho_grp}->{'species'}->{$o_species}->{'gene'} };
        my @o_genes      = ();
        foreach my $o_gene (@orig_o_genes) {
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

