#!/usr/bin/env perl

# build_func_boolean.pl -- Erich Schwarz <ems394@cornell.edu>, 3/2/2016.
# Purpose: build FUNC-usable Boolean from file of positive genes (one per line) and file of all genes to be considered; read gene names at start of each line, ignoring rest of text; reject redundant names, or positives missing from general list; ignore comment lines starting with '#'.  Optionally, tolerate lists of genes with non-standard gene names (but do not allow them into the output pre-FUNC file).

use strict;
use warnings;
use Getopt::Long;

my $positive_genes = q{};
my $general_genes  = q{};
my $help;

my $gene = q{};
my $data_ref;
my $tolerate;

GetOptions ( 'positive=s' => \$positive_genes,
             'general=s'  => \$general_genes,
             'tolerate'   => \$tolerate,
             'help'       => \$help,           );

if ($help or (! $positive_genes) or (! $general_genes) ) { 
    die "Format: build_func_boolean.pl\n",
        "            -p|--positive  [genes to get '1']\n",
        "            -g|--general   [*all* genes, defaulting to '0']\n",
        "            -t|--tolerate  [loudly reject non-standard gene names in subset gene lists, but do not die]\n",
        "            -h|--help\n",
        "        [print to STDOUT]\n",
        ;
}

open my $GENERAL, '<', $general_genes or die "Can't open file with general gene list $general_genes: $!";
while (my $input = <$GENERAL>) { 
    chomp $input;
    # Silently pass over comment lines.
    if ( $input !~ /\A \# /xms ) { 
        if ( $input =~ /\A (\S+) /xms ) { 
            $gene = $1;
            if ( exists $data_ref->{'gene'}->{$gene} ) { 
                die "Redundant gene $gene in general gene list!\n";
            }
            $data_ref->{'gene'}->{$gene} = 1;
        }
    }
}
close $GENERAL or die "Can't close filehandle to general gene list: $!";

open my $POSITIVE, '<', $positive_genes or die "Can't open file with positive gene list $positive_genes: $!";
while (my $input = <$POSITIVE>) {
    chomp $input;
    # Silently pass over comment lines.
    if ( $input !~ /\A \# /xms ) {
        if ( $input =~ /\A (\S+) /xms ) {
            $gene = $1;
            if (! exists $data_ref->{'gene'}->{$gene} ) {
                warn "Gene $gene in positive but not general gene list!\n";
                if (! $tolerate) {
                    die "\n";
                }
            }
            if ( exists $data_ref->{'positive'}->{$gene} ) { 
                die "Redundant gene $gene in positive gene list!\n";
            }
            # Nonstandard names can be tolerated, but not allowed in downstream FUNC lists.
            if ( exists $data_ref->{'gene'}->{$gene} ) {
                $data_ref->{'positive'}->{$gene} = 1; 
            }
        }
    }
}
close $POSITIVE or die "Can't close filehandle to general gene list: $!";

foreach my $gene1 (sort keys %{ $data_ref->{'gene'} }) { 
    if ( exists $data_ref->{'positive'}->{$gene1} ) { 
        print "$gene1\t1\n";
    }
    else { 
        print "$gene1\t0\n";
    }
}

