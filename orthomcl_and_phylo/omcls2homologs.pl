#!/usr/bin/env perl

# omcls2homologs.pl -- Erich Schwarz <emsch@caltech.edu>, 12/22/2012.
# Purpose: given OrthoMCL report text (file or piped), a primary taxon, and a secondary taxon, annotate the primary taxon's genes/proteins with either strict or any orthologs from the secondry taxon.

use strict;
use warnings;
use Getopt::Long;

my @infiles         = ();
my $primary_taxon   = q{};
my $secondary_taxon = q{};
my $require_strict;
my $help;

my $data_ref;

GetOptions ( 'infile=s{,}'    => \@infiles,
             'primary=s'      => \$primary_taxon,
             'secondary=s'    => \$secondary_taxon,
             'require_strict' => \$require_strict,
             'help'           => \$help,    );

if ( $help or (! $primary_taxon ) or (! $secondary_taxon) ) {
    die "Format: omcls2homologs.pl\n",
        "        --primary|-p         [primary taxon, whose genes/proteins are the index for orthologies]\n",
        "        --secondary|-s       [secondary taxon, whose orthologous genes/proteins are tabulated either as strict or as any orthologs, by name]\n",
        "        --infile|-i          [input files or '-' for stream]\n",
        "        --require_strict|-r  [only print lines where a strict ortholog exists]\n",
        "        --help|-h            [print this message]\n",
        ;
}

# To be used as filehandle for the input(s):
my $INFILE;

foreach my $infile (@infiles) { 
    if ($infile eq '-') {
        # Special case: get the stdin handle
        $INFILE = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $INFILE, '<', $infile or die "Can't open input file $infile. $!\n";
    }
    while (my $input = <$INFILE>) { 
        chomp $input;

        # Sample input line:  
        # ORTHOMCL19007(2 genes,2 taxa: WBGene00031156(briggsae) WBGene00005532(elegans)

        if (     ( $input =~ /\S/xms ) 
             and ( $input !~ / \A 
                               ORTHOMCL\d+ 
                               \( \d+ \s+ genes,\d+ \s+ taxa\) : 
                               (?: \s+ [^\s\(\)]+ \( [^\s\(\)]+ \) )+
                               \s* \z /xms ) 
           ) { 
            die "Obviously not an OrthoMCL text line!\n";
        }
        else { 
            my @primary_genes   = ();
            my @secondary_genes = ();
            while ( $input =~ / (\S+) \( $primary_taxon \) /gxms ) { 
                my $primary_gene = $1;
                push @primary_genes, $primary_gene;
            }
            while ( $input =~ / (\S+) \( $secondary_taxon \) /gxms ) {
                my $secondary_gene = $1;
                push @secondary_genes, $secondary_gene;
            }
            foreach my $primary_gene1 (@primary_genes) { 
                foreach my $secondary_gene1 (@secondary_genes) { 
                    if ( exists $data_ref->{'primary_gene'}->{$primary_gene1}->{'secondary_gene'}->{$secondary_gene1} ) { 
                        warn "Predicted orthology between primary gene $primary_gene1 and secondary gene $secondary_gene1 observed in more than one instance.\n";
                    }
                    # Record *both* directions, so that I can later ensure that 'strict' means 1:1, and not 2:1 or 1:2.
                    $data_ref->{'primary_gene'}->{$primary_gene1}->{'secondary_gene'}->{$secondary_gene1} = 1;
                    $data_ref->{'secondary_gene'}->{$secondary_gene1}->{'primary_gene'}->{$primary_gene1} = 1;
                }
            }
        }
    }
    close $INFILE or die "Can't close filehandle to input file $infile. $!\n";
}

my $header = "Gene\tStrict_orthologs\tOrthologs\n";

my @primary_genes = sort keys %{ $data_ref->{'primary_gene'} };

foreach my $primary_gene (@primary_genes) { 
    my @secondary_genes                  = ();
    my $secondary_ortholog_count         = 0;
    my $strict_ortholog                  = q{};
    my $ortholog_list                    = q{};

    if ( exists $data_ref->{'primary_gene'}->{$primary_gene}->{'secondary_gene'} ) { 
        @secondary_genes = sort keys %{ $data_ref->{'primary_gene'}->{$primary_gene}->{'secondary_gene'} };
    }
    $secondary_ortholog_count = @secondary_genes;

    if ( $secondary_ortholog_count == 1 ) { 
        my @recursive_primary_orthologs = sort keys %{ $data_ref->{'secondary_gene'}->{$secondary_genes[0]}->{'primary_gene'} };
        my $recursive_primary_ortholog_count = @recursive_primary_orthologs;
        if ( $recursive_primary_ortholog_count == 1 ) {
            $strict_ortholog = $secondary_genes[0];
        }
    }

    if ( $secondary_ortholog_count >= 1 ) { 
        $ortholog_list = join '; ', @secondary_genes;
    }

    print $header if $header;
    $header = q{};

    if ( $strict_ortholog or (! $require_strict) ) { 
        print "$primary_gene\t$strict_ortholog\t$ortholog_list\n";
    }
}

