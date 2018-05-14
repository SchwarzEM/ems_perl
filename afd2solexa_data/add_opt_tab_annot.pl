#!/usr/bin/env perl

# add_opt_tabannot.pl -- Erich Schwarz <emsch@caltech.edu>, 7/10/2012.
# Purpose: given file with gene names, list to annot, and annot, print file w/ tab +/- annot.

use strict;
use warnings;
use Getopt::Long;

my $data_file = q{};
my $genelist  = q{};
my $annot     = q{};
my %gene2annot = ();

my $help;

GetOptions ( 'data=s'  => \$data_file,
             'genes=s' => \$genelist,
             'annot=s' => \$annot,
             'help'    => \$help, );

if ( $help or (! $data_file ) or (! $genelist ) or (! $annot ) ) { 
    die "Format: add_opt_tab_annot.pl\n",
        "            --data|-d  <data_file to recolumn>\n",
        "            --genes|-g <gene_list>\n",
        "            --annot|-a <text to annotate to select genes>\n",
        "            --help|-h \n",
        ;
}

open my $GENES, '<', $genelist or die "Can't open gene list $genelist: $!";
while (my $gene = <$GENES>) { 
    chomp $gene;
    if ( $gene =~ /\S/xms ) { 
        $gene2annot{$gene} = 1;
    }
}
close $GENES or die "Can't close filehandle to gene list $genelist: $!";

open my $DATA, '<', $data_file or die "Can't open data file $data_file: $!";
while (my $input = <$DATA>) { 
    chomp $input;
    $input = $input . "\t";
    if ( $input =~ /\A (\S+) /xms ) { 
        my $gene = $1;
        if ($gene2annot{$gene}) { 
            $input = $input . $annot;
        }
    }
    print "$input\n";
}
close $DATA or die "Can't close filehandle to data file $data_file: $!";

