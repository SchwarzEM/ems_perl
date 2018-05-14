#!/usr/bin/env perl

# xatab2wbannot.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/18/2010.
# Purpose: given TableMaker outputs of various sorts for WBGenes, export one-line-per-gene annotations.

use strict;
use warnings;
use Getopt::Long;

my $data   = q{};
my $gene   = q{};
my $annot  = q{};
my $annot1 = q{};
my $annot2 = q{};

my $gene_info_ref;

my %OK_data = ( KOG  => 1,
                desc => 1, );

my $help;

GetOptions ( 'data=s', => \$data,
             'help'    => \$help, );

if ($help or (! $data) ) { 
    die "Format: xatab2wbannot.pl",
        "    --data|-d [KOG|desc]",
        " <stream or input files>",
        " --help|-h\n",
        ;
}

if (! $OK_data{$data}) { 
    die "The data type \"$data\" is not supported.\n";
}

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ / \A 
                     \" (WBGene\d+) \" 
                     .* 
                     \t \" ([^\t\"]+) \" 
                     \t \" ([^\t\"]+) \"
                     \z /xms ) { 
        if ($data eq 'KOG') { 
            $gene   = $1;
            $annot1 = $2;
            $annot2 = $3;
            $annot  = "$annot1 ($annot2)";
            $gene_info_ref->{$gene}->{$annot} = 1;
        }
        if ($data eq 'desc') { 
            $gene   = $1;
            $annot  = $3;
            # Get rid of some unavoidable ACeDB text cruft:
            $annot =~ s/\\//g;
            $gene_info_ref->{$gene}->{$annot} = 1;
        }
    }
}

foreach my $wb_gene (sort keys %{ $gene_info_ref }) { 
    my @annots = ();
    @annots = sort keys %{ $gene_info_ref->{$wb_gene} };
    $annot = join "; ", @annots;
    print "$wb_gene\t\"$annot\"\n";
}

