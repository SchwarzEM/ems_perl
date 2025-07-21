#!/usr/bin/env perl

# get_largest_isoforms.pl -- Erich Schwarz <ems394@cornell.edu>, 7/29/2018.
# Purpose: given a proteome in FASTA format with known pattern of identifiers in the headers, print largest isoform per gene (currently handles WormBase official, augustus, augustus-like '.t\d+', ensembl, FlyBase, Maker (older and newer formats), ParaSite [release 10], ParaSite [pre-release 10], or three-column proteomes).

use strict;
use warnings;
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

                   maker_old => 'maker_old',
                   mak_old   => 'maker_old',

                   parasite => 'parasite',
                   par      => 'parasite',

                   parasite_old => 'parasite_old',
                   par_old      => 'parasite_old',

                   liftoff  => 'liftoff',

                   column3  => 'column3',
                   col3     => 'column3',

                   irregular1 => 'irregular1',
                   irr1       => 'irregular1',
);

my $header  = q{};
my $protein = q{};
my $gene    = q{};

my $help;

GetOptions ( 'infiles=s{,}' => \@infiles,
             'type=s'       => \$header_type,
             'help'         => \$help,   );

if ( $help or (! @infiles) ) { 
    die "Format: get_largest_isoforms.pl\n",
        "    --infile|-i   <input stream/files>\n",
        "    --type|-t     [header types:\n",
        "                   'wormbase|wb'\n",
        "                   'augustus|aug'\n",
        "                   'aug_like'\n",
        "                   'ensembl|ens'\n",
        "                   'flybase|fly'\n",
        "                   'fly_old'\n",
        "                   'maker|mak'\n",
        "                   'maker_old|mak_old'\n",
        "                   'parasite|par'\n",
        "                   'parasite_old|par_old'\n",
        "                   'liftoff'\n",
        "                   'column3|col3';\n", 
        "                or 'irregular1|irr3';\n",
        "                   default is 'wormbase']\n",
        "    --help|-h     [print this message]\n",
        ;
}

if ($header_type) { 
    if (! exists $ok_headers{$header_type} ) { 
        die "Header type must be",
            " 'wormbase|wb', 'augustus|aug', 'aug_like', 'ensembl|ens', 'flybase|fly',",
            " 'fly_old', 'maker|mak', 'maker_old|mak_old', 'parasite|par', 'parasite_old|par_old',",
            " or 'column3|col3', not \"$header_type\"\n",
            ;
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

            if ( ( $header_type eq 'wormbase' ) and ( $input =~ /\A > ( (\S+) \s .* (WBGene\d+) .*) \z/xms ) ) { 
                $header  = $1;
                $protein = $2;
                $gene    = $3;
                $data_ref->{'gene'}->{$gene}->{'protein'}->{$protein}->{'header'} = $header;
            }

            elsif ( ( $header_type eq 'ensembl' ) and ( $input =~ /\A > ( (\S+) \b .* \s gene[:] (\S+) \b .* ) \z/xms ) ) {
                $header  = $1;
                $protein = $2;
                $gene    = $3;
                $data_ref->{'gene'}->{$gene}->{'protein'}->{$protein}->{'header'} = $header;
            }

            elsif ( ( $header_type eq 'parasite' ) and ( $input =~ /\A > ( (\S+) \b .* \s gene= (\S+) \b .* ) \z/xms ) ) {
                $header  = $1;
                $protein = $2;
                $gene    = $3;
                $data_ref->{'gene'}->{$gene}->{'protein'}->{$protein}->{'header'} = $header;
            }      

            elsif ( ( $header_type eq 'parasite_old' ) and ( $input =~ /\A > ( (\S+) \b .* \s gene_id= (\S+) \b .* ) \z/xms ) ) {
                $header  = $1;
                $protein = $2;
                $gene    = $3;
                $data_ref->{'gene'}->{$gene}->{'protein'}->{$protein}->{'header'} = $header;
            }

            elsif ( ( $header_type eq 'liftoff' ) and ( $input =~ /\A > ( (\S+) \b .* gene_id= ([^;\s]+) [;] .* ) \z/xms ) ) {
                $header  = $1;
                $protein = $2;
                $gene    = $3;
                $data_ref->{'gene'}->{$gene}->{'protein'}->{$protein}->{'header'} = $header;
            }


            elsif ( ( $header_type eq 'augustus' ) and ( $input =~ /\A > ( ( (\S+ \.g\d+) \.t\d+ ) \b .* ) \z/xms ) ) {
                $header  = $1;
                $protein = $2;
                $gene    = $3;
                $data_ref->{'gene'}->{$gene}->{'protein'}->{$protein}->{'header'} = $header;
            }

            elsif ( ( $header_type eq 'augustus_like') and ( $input =~ /\A > ( ( (\S+) \.t\d+ ) \b .* ) \z/xms ) ) {
                $header  = $1;
                $protein = $2;
                $gene    = $3;
                $data_ref->{'gene'}->{$gene}->{'protein'}->{$protein}->{'header'} = $header;
            }

            elsif ( ( $header_type eq 'flybase' ) and ( $input =~ /\A > ( (\S+) .* parent=(FBgn\d+) .+ ) \z/xms ) ) {
                $header  = $1;
                $protein = $2;
                $gene    = $3;
                $data_ref->{'gene'}->{$gene}->{'protein'}->{$protein}->{'header'} = $header;
            }

            elsif ( ( $header_type eq 'flybase_old' ) and ( $input =~ /\A > ( ( (\S+) [-]R[A-Z] ) \b .*) \z/xms ) ) {
                $header  = $1;
                $protein = $2;
                $gene    = $3;
                $data_ref->{'gene'}->{$gene}->{'protein'}->{$protein}->{'header'} = $header;
            }

            elsif ( ( $header_type eq 'maker' ) and ( $input =~ /\A > ( ( (\S+) [-]\d+) \b .*) \z/xms ) ) {
                $header  = $1;
                $protein = $2;
                $gene    = $3; 
                $data_ref->{'gene'}->{$gene}->{'protein'}->{$protein}->{'header'} = $header;
            }

            elsif ( ( $header_type eq 'maker_old' ) and ( $input =~ /\A > ( ( (\S+) [-]mRNA[-]\d+) \b .*) \z/xms ) ) {
                $header  = $1;
                $protein = $2;
                $gene    = $3; 
                $data_ref->{'gene'}->{$gene}->{'protein'}->{$protein}->{'header'} = $header;
            }

            elsif ( ( $header_type eq 'column3' ) and ( $input =~ /\A > ( (\S+) \s+ \S+ \s+ (\S+) \s* ) \z/xms ) ) {
                $header  = $1;
                $protein = $2;
                $gene    = $3;
                $data_ref->{'gene'}->{$gene}->{'protein'}->{$protein}->{'header'} = $header;
            }

            elsif ( ( $header_type eq 'irregular1' ) and ( $input =~ /\A > ( (\S+) .*) \z /xms ) ) {
                $header  = $1;
                $protein = $2;
                # Default is for each protein name to exactly match its gene name:
                $gene    = $protein;
                # Alternatively, if the protein has a ".\d+" suffix:
                if ( $protein =~ /\A (\S+) \. \d+ \z/xms ) {
                    $gene = $1;
                }
                $data_ref->{'gene'}->{$gene}->{'protein'}->{$protein}->{'header'} = $header;
            }

            else { 
                die "Can't parse header: $input\n";
            }
        }
        else { 
            if ( ( $input =~ /\S/xms ) and ( $gene and $protein ) ) { 
                $data_ref->{'gene'}->{$gene}->{'protein'}->{$protein}->{'sequence'} .= $input;
            }
        }
    }
    close $INPUT_FILE or die "Can't close filehandle to input file $infile. $!\n";
}

# List all the genes observed in the proteome:
my @genes = sort keys %{ $data_ref->{'gene'} };

# Print the longest isoform for each gene in turn:
foreach my $gene1 (@genes) { 
    # Select its largest isoform, with a background default sort of ASCII-betical text:
    my @isoforms = sort {     length( $data_ref->{'gene'}->{$gene1}->{'protein'}->{$b}->{'sequence'} ) 
                          <=> length( $data_ref->{'gene'}->{$gene1}->{'protein'}->{$a}->{'sequence'} ) } 
                   sort keys %{ $data_ref->{'gene'}->{$gene1}->{'protein'} };
    my $largest_isoform = $isoforms[0];

    # Print its header:
    print '>', $data_ref->{'gene'}->{$gene1}->{'protein'}->{$largest_isoform}->{'header'}, "\n";

    # Print its FASTA-formatted sequence text:
    my $curr_seq = $data_ref->{'gene'}->{$gene1}->{'protein'}->{$largest_isoform}->{'sequence'};
    my @output_lines = unpack( "a60" x ( length($curr_seq)/60 + 1), $curr_seq );
    foreach my $output_line (@output_lines) {
        print "$output_line\n" if ($output_line =~ /\S/);
    }
}

