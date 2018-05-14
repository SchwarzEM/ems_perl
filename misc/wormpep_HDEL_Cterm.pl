#!/usr/bin/env perl

# wormpep_HDEL_Cterm.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/14/2009.
# Purpose: scan wormpepXX for proteins ending in HDEL.

use strict;
use warnings;

my $cds  = q{};
my $gene = q{};
my $seq  = q{};

while (my $input = <>) { 
    chomp $input;
    # New FASTA sequence header:
    if ( $input =~ / \A > (\S+) \s* /xms ) { 
            # Capture new value:
            my $new_cds = $1;

            # Process stored older values:
            evaluate_data_and_print_results($cds,$gene,$seq);

            # Then reset values and start recording data again:
            $seq  = q{};
            $cds = $new_cds;
            $gene = $cds;
            $gene =~ s/[a-zA-Z]\z//;
    }
    # Three-letter gene name:
    if ( $input =~ / \A > \S+ \s+ .* locus: (\S+) \s /xms ) { 
        $gene = $1;
    }
    # Plain amino acid sequence data:
    if ( $input =~ / \A [a-zA-Z]+ \z /xms ) { 
        $seq .= $input;
    }
}
evaluate_data_and_print_results($cds,$gene,$seq);

sub evaluate_data_and_print_results { 
    my $_cds  = $_[0];
    my $_gene = $_[1];
    my $_seq  = $_[2];
    if ( $_seq =~ / HDEL \z /ixms ) { 
        print "$_gene\t$_cds\n";
    }
}


