#!/usr/bin/env perl

# strict_vs_exp_omcls_legacy.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/16/2009.
# Purpose: [LEGACY:] given an OrthoMCL line, allow one or more taxa to expand; require others be 0-1 genes.

use strict;
use warnings;
use Getopt::Long;

my $input_omcl = q{};

my @variable   = ();
my %var_taxa   = ();

my @expandable = ();
my %exp_taxa   = ();

my $help;

GetOptions ( 'input_omcl:s'     => \$input_omcl,
             'variable=s{1,}'   => \@variable,
             'expandable=s{1,}' => \@expandable,
             'help'             => \$help,       );

$input_omcl ||= 'all_orthomcl.out';   # Default name, if not specified.

if ( ($help) or (! -r $input_omcl ) ) { 
     die "Format: strict_vs_exp_omcls.pl --input_omcl|-i [orthomcl.out] --variable|-v [allow 0-N genes for these taxa] --expandable|-e [allow 2+ genes for these taxa]\n";
}

foreach my $ok_taxon (@expandable) { 
    $exp_taxa{$ok_taxon} = 1;
}

foreach my $var_taxon (@variable) {
    $var_taxa{$var_taxon} = 1;
}

open my $INPUT_OMCL, '<', $input_omcl or die "Cannot open OrthoMCL file $input_omcl: $!";
while (my $input = <$INPUT_OMCL>) { 
    chomp $input;
    my $printable = 1;
    my %seen_taxon = ();
    # Enforce correct format:
    if ($input !~ / \A
                    ORTHOMCL\d+
                    \( \d+\sgenes,\d+\staxa \)
                    : \s+
                    .+\S
                    \s*  
                    \z     
                  /xms) { 
        die "From OrthoMCL file $input_omcl, can't parse input line: $input\n";
    }

    if ($input =~ / \A
                    ORTHOMCL\d+ 
                    \( \d+\sgenes,\d+\staxa \)    # just (punctuation)
                    : \s+ 
                    (.+\S)                        # $1 -> $oprots_line
                    \s* 
                    \z 
                  /xms) { 
        my $oprots_line = $1;
        my @orthoprots = split /\s+/, $oprots_line;
        foreach my $o_prot (@orthoprots) { 
            # Enforce correct format:
            if ( $o_prot !~ / \A [^\s\(\)]+ \( [^\s\(\)]+ \) \z /xms ) { 
                die "Can't parse protein $o_prot in input $input!\n";
            }

            # If format's OK, extract, check and store various mappings:
            if ( $o_prot =~ / \A [^\s\(\)]+ \( ( [^\s\(\)]+ ) \) \z /xms ) { 
                my $species = $1;
                # Only allow an unflagged species to be seen once:
                if ( ( (! exists $exp_taxa{$species} ) and (! exists $var_taxa{$species} ) ) and ( exists $seen_taxon{$species} ) ) { 
                    $printable = 0;
                }
                $seen_taxon{$species} += 1;
            }
        }
    }

    # Require that expanded taxa have >=2 representatives in the orthology group:
    foreach my $exp_taxon (sort keys %exp_taxa) { 
        if ( (! $seen_taxon{$exp_taxon} ) or ( $seen_taxon{$exp_taxon} < 2 ) ) {
            $printable = 0;
        }     
    }

    # Finally, print any lines that passed the above tests.
    if ($printable) { 
        print "$input\n";
    }
}

