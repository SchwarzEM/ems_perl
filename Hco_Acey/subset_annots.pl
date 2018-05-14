#!/usr/bin/env perl

# subset_annots.pl -- Erich Schwarz <emsch@caltech.edu>, 11/16/2012.
# Purpose: given a list of genes and a file annotating them (along with others), print the subset of the file with the listed genes.  First devised to cope with RSEM slices for Acey.

use strict;
use warnings;
use Getopt::Long;

my $genelist = q{};
my $annots   = q{};

my $data_ref;

my $help;

GetOptions ( 'genelist=s' => \$genelist,
             'annots=s'   => \$annots,
             'help'       => \$help, );

if ( $help or (! $genelist) or (! $annots) ) { 
    die "Format: subset_annots.pl --genelist|-g [list of genes to select] --annots|-a [annotation file to select from] --help|-h\n";
}

open my $GENES, '<', $genelist or die "Can't open gene list: $genelist\n";
while (my $input = <$GENES>) { 
    chomp $input;
    $data_ref->{'ok_gene'}->{$input} = 1;
}
close $GENES or die "Can't close filehandle to gene list: $genelist\n";

open my $ANNOTS, '<', $annots or die "Can't open annotations file: $annots\n";
while (my $input = <$ANNOTS>) { 
    chomp $input; 
    if ( $input =~ /\A (\S+) \s+ /xms ) { 
        my $gene = $1;
        if ( exists $data_ref->{'ok_gene'}->{$gene} ) { 
            print "$input\n";
         }
    }
}
close $ANNOTS or die "Can't close filehandle to annotations file: $annots\n";

