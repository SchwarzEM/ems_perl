#!/usr/bin/perl

# extract_psi-blast_hitlist.pl -- Erich Schwarz, 5/5/2008 (orig. 6/28/02).
# Purpose: get hitlists from last round of psi-Blast (convergence not req.).

use strict;
use warnings;

# Require a single input file to keep data clean:

my $input_file = $ARGV[0];

if ( (! $input_file) or ($#ARGV > 0) ) { 
    die "Format: ./extract_psi-blast_hitlist.pl [single input file].\n";
}

open my $INFILE, '<', $input_file 
    or die "Can't open input file $input_file: $!";

# [Format of things to extract:]

# Sequences producing significant alignments:                      (bits) Value
# Sequences used in model and found again:

# E02H4.1 CE05547  locus:del-1 mechanosensory protein status:Confi...   397   e-111
# T01C8.7 CE25094  locus:mec-4  status:Partially_confirmed              379   e-105
# C47C12.6 CE24854  locus:deg-1  status:Confirmed                       376   e-105

# Sequences not found previously or not previously below threshold:

# F59F3.4 CE11544    status:Predicted TR:Q93808 protein_id:CAA91993.1    37   0.007
# C47F8.7 CE17576    status:Predicted TR:O62110 protein_id:CAA15834.1    29   2.5  
# F26A3.6 CE09671    status:Predicted TR:Q93597 protein_id:CAB01704.1    28   4.5  
# AC3.3 CE05133  locus:pqn-1  status:Predicted TR:Q17400 protein_i...    28   7.0  
# B0334.11 CE19661  locus:ooc-3  status:Partially_confirmed TR:Q9X...    27   8.9  
# 
# CONVERGED!

my %list_of_hits = ();
my $scan_hits = 0;

while (my $input = <$INFILE>) {
    chomp $input;

    # If a new round of psi-Blast hits is starting:
    if ( $input =~ / \A Sequences[ ] 
                        used[ ] 
                        in[ ] 
                        model[ ] 
                        and[ ] 
                        found[ ] 
                        again: /xms ) {

        # Zero out the hitlist.
        %list_of_hits = ();

        # And start reading new seqs.
        $scan_hits = 1;
    }

    # When list of confirmed positives stops...
    elsif ( $input =~ / \A 
                        Sequences[ ] 
                        not[ ] 
                        found[ ] 
                        previously[ ] 
                        or[ ] 
                        not[ ] 
                        previously[ ] 
                        below[ ] 
                        threshold: /xms ) {
        # ... stop reading in seqs.
        $scan_hits = 0;
    }

    # N.B.: this script has no specific requirement for convergence.

    # But if reading hits, well, just read 'em:
    elsif ( ( $scan_hits ) 
             and ( $input =~ / \A (\S+) \s* /xms ) ) { 
        my $seq = $1;
        $list_of_hits{$seq} = 1;
    }
}

close $INFILE 
    or die "Can't close filehandle of input file $input_file: $!";

foreach my $seq (sort keys %list_of_hits) {
        print "$seq\n";
}

