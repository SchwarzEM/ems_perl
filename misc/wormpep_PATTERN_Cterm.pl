#!/usr/bin/env perl

# wormpep_PATTERN_Cterm.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/14/2009.
# Purpose: scan wormpepXX for proteins ending in PATTERN.

use strict;
use warnings;
use Getopt::Long;

my $cds     = q{};
my $gene    = q{};
my $seq     = q{};
my $PATTERN = q{};

GetOptions ( "pattern=s" => \$PATTERN, );

if ( (! $PATTERN ) 
     or ( $PATTERN !~ / \A \S+ \z /xms ) ) { 
    die 'Format: ./wormpep_PATTERN_Cterm.pl --pattern|-p [\S+] <input files or STDIN>', "\n";
}

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
    if ( $_seq =~ / $PATTERN \z /xmis ) { 
        print "$_gene\t$_cds\n";
    }
}

