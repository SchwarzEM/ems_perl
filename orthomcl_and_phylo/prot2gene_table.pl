#!/usr/bin/env perl

# prot2gene_table.pl -- Erich Schwarz <ems394@cornell.edu>, 3/21/2020.
# Purpose: given a set of proteome files, tabulate CDSes-to-genes, with none in >=2 proteomes, and with optional third column with a taxon name; prerequisite for prot2gene_omcl.pl.

# Note 1: this uses syntax taken straight form get_largest_isoforms.pl -- Erich Schwarz <ems394@cornell.edu>, 7/29/2018. 
# Note 2: this maps CDSes to unique genes; if it mapped *proteins*, the mapping could not be unique
#     because it is entirely possible for a single polypeptide to be encoded by >2 genes.

# Import proteome files by '--argument [file list]'.
#     Each argument groups proteomes; each group is then assigned 
#     a parsing rule.  This allows modular addition of new parsing 
#     rules, as needed ad-hoc, within a stable frame of general, 
#     unchanging code.

use strict;
use warnings;
use autodie;
use Getopt::Long;

my $data_ref;

my @infiles     = ();
my $header_type = q{};

my %ok_headers = ( wormbase => 'wormbase',
                   wb       => 'wormbase',

                   augustus => 'augustus',
                   aug      => 'augustus',

                   aug_like => 'augustus_like',

                   ensembl  => 'ensembl',
                   ens      => 'ensembl',

                   flybase  => 'flybase',
                   fly      => 'flybase',

                   fly_old  => 'flybase_old',

                   maker    => 'maker',
                   mak      => 'maker', 

                   parasite => 'parasite',
                   par      => 'parasite',

                   par_old  => 'parasite_old',

                   column3  => 'column3',
                   col3     => 'column3',
);

my $header  = q{};
my $protein = q{};
my $gene    = q{};
my $species = q{};

my $help;

GetOptions ( 'infiles=s{,}' => \@infiles,
             'type=s'       => \$header_type,
             'species=s'    => \$species,
             'help'         => \$help,   );

if ( $help or (! @infiles) ) { 
    die "Format: prot2gene_table.pl\n",
        "    --infile|-i   <input stream/files>\n",
        "    --type|-t     [header types:\n",
        "                   'wormbase|wb'\n",
        "                   'augustus|aug'\n",
        "                   'aug_like'\n",
        "                   'ensembl|ens'\n",
        "                   'flybase|fly'\n",
        "                   'fly_old'\n",
        "                   'maker|mak'\n",
        "                   'parasite|par'\n",
        "                   'par_old'\n",
        "                or 'column3|col3';\n", 
        "                   default is 'wormbase']\n",
        "    --species|-s  [optional: provide a taxon name, to be printed as a third column]\n",
        "    --help|-h     [print this message]\n",
        ;
}

if ($header_type) { 
    if (! exists $ok_headers{$header_type} ) { 
        die "Header type must be 'wormbase|wb', 'augustus|aug', 'aug_like', 'ensembl|ens', 'flybase|fly', 'fly_old', 'maker|mak', 'parasite|par', 'par_old', or 'column3|col3', not \"$header_type\"\n";
    }
    else { 
        # Map to full names of header types, if abbreviated:
        $header_type = $ok_headers{$header_type};
    }
}

foreach my $infile (@infiles) { 
    my $INPUT_FILE;
    if ($infile eq '-') {
        # Special case: get the stdin handle
        $INPUT_FILE = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $INPUT_FILE, '<', $infile or die "Can't open input file $infile. $!\n";
    }
    while (my $input = <$INPUT_FILE>) { 
        chomp $input;
        if ( $input =~ /\A > /xms ) { 

            # As long as protein-gene relations can be reliably extracted from headers with regexes, we can easily add all the rules we want.
            # We just need to get $header, $protein, and $gene reliably from the header lines *somehow*.
            # So this part can be expanded easily; after this step, the rest of the script is changeless.

            if ( ( $header_type eq 'wormbase' ) and ( $input =~ /\A > (\S+) \s .* (WBGene\d+) .* \z/xms ) ) { 
                $protein = $1;
                $gene    = $2;
                $data_ref->{'protein'}->{$protein}->{'gene'} = $gene;
            }

            elsif ( ( $header_type eq 'ensembl' ) and ( $input =~ /\A > (\S+) \b .* \s gene[:] (\S+) \b .* \z/xms ) ) {
                $protein = $1;
                $gene    = $2;
                $data_ref->{'protein'}->{$protein}->{'gene'} = $gene;
            }

            elsif ( ( $header_type eq 'parasite' ) and ( $input =~ /\A > (\S+) \b .* \s gene= (\S+) \b .* \z/xms ) ) {
                $protein = $1;
                $gene    = $2;
                $data_ref->{'protein'}->{$protein}->{'gene'} = $gene;
            }

            elsif ( ( $header_type eq 'parasite_old' ) and ( $input =~ /\A > (\S+) \b .* \s gene_id= (\S+) \b .* \z/xms ) ) {
                $protein = $1;
                $gene    = $2;
                $data_ref->{'protein'}->{$protein}->{'gene'} = $gene;
            }

            elsif ( ( $header_type eq 'augustus' ) and ( $input =~ /\A > ( (\S+ \.g\d+) \.t\d+ ) \b .* \z/xms ) ) {
                $protein = $1;
                $gene    = $2;
                $data_ref->{'protein'}->{$protein}->{'gene'} = $gene;
            }

            elsif ( ( $header_type eq 'augustus_like') and ( $input =~ /\A > ( (\S+) \.t\d+ ) \b .* \z/xms ) ) {
                $protein = $1;
                $gene    = $2;
                $data_ref->{'protein'}->{$protein}->{'gene'} = $gene;
            }

            elsif ( ( $header_type eq 'flybase' ) and ( $input =~ /\A > (\S+) .* parent=(FBgn\d+) .+ \z/xms ) ) {
                $protein = $1;
                $gene    = $2;
                $data_ref->{'protein'}->{$protein}->{'gene'} = $gene;
            }

            elsif ( ( $header_type eq 'flybase_old' ) and ( $input =~ /\A > ( (\S+) [-]R[A-Z] ) \b .* \z/xms ) ) {
                $protein = $1;
                $gene    = $2;
                $data_ref->{'protein'}->{$protein}->{'gene'} = $gene;
            }

            elsif ( ( $header_type eq 'maker' ) and ( $input =~ /\A > ( (\S+) [-]mRNA[-]\d+) \b .* \z/xms ) ) {
                $protein = $1;
                $gene    = $2;
                $data_ref->{'protein'}->{$protein}->{'gene'} = $gene;
            }

            elsif ( ( $header_type eq 'column3' ) and ( $input =~ /\A > (\S+) \s+ \S+ \s+ (\S+) \s* \z/xms ) ) {
                $protein = $1;
                $gene    = $2;
                $data_ref->{'protein'}->{$protein}->{'gene'} = $gene;
            }

            else { 
                die "Can't parse header: $input\n";
            }
        }
    }
    close $INPUT_FILE or die "Can't close filehandle to input file $infile. $!\n";
}

# List all the proteins observed in the proteome:
my @proteins = sort keys %{ $data_ref->{'protein'} };

# Print a table with each protein, its gene, and optionally its species/taxon:
foreach my $protein1 (@proteins) {
    my $gene1 = $data_ref->{'protein'}->{$protein1}->{'gene'} ;
    print "$protein1\t$gene1";
    if ( $species ) {
        print "\t$species";
    }
    print "\n";
}

