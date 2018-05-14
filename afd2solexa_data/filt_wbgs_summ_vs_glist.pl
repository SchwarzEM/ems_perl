#!/usr/bin/perl

# filt_wbgs_summ_vs_glist.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/22/2008.
# Purpose: given a plain genelist, either accept or reject wbgenes_summ genes that are found in the genelist.

use strict;
use warnings;

# Enforce correct argument count and first argument.

if ($#ARGV != 2) { 
    &fail_loudly;
} 

my $accept = 0;
my $direction = shift @ARGV;
unless ( ($direction eq '--accept') or ($direction eq '--reject') ) { 
    &fail_loudly;
}
if ($direction eq '--accept') { 
    $accept = 1;
}
if ($direction eq '--reject') { 
    $accept = 0;
}

sub fail_loudly { 
    die 'Format: ./filt_wbgs_summ_vs_glist.pl (--accept|--reject) ',
        '[genelist] [wbgenes_summ.txt]',
        "\n",
    ;
}

# Read in genelist.

my %listed_genes = ();

my $genelist = $ARGV[0];

open my $GENES, "<", "$genelist" 
    or die "Can't open gene list $genelist\n";

while (my $input = <$GENES>) { 
    chomp $input;
    if ($input =~ / \A ( WBGene\d+\|\S+ ) \s* \z /xms) {
        $listed_genes{$1} = 1;
    }
}

close $GENES;

# Filter the wbgenes_summ genes.

my $wb_summ = $ARGV[1];

open my $WB_SUMM, "<", "$wb_summ" 
    or die "Can't open wbgenes_summ file $wb_summ\n";

while (my $input = <$WB_SUMM>) { 
    chomp $input;
    my $wbs_gene = q{};
    if ($input =~ / \A ( WBGene\d+\|\S+ ) \s* /xms) {   # was '\s'
        $wbs_gene = $1;
    }
    if ( $wbs_gene ) { 
        if ( ( int($accept) == 1 ) and (  $listed_genes{$wbs_gene} ) ) {
            print "$input\n";
        }
        if ( ( int($accept) == 0 ) and (! $listed_genes{$wbs_gene} ) ) {
            print "$input\n";
        }
    }
}

