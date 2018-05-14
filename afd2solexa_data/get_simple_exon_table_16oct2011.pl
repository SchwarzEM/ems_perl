#!/usr/bin/env perl

# get_simple_exon_table_16oct2011.pl -- Erich Schwarz <emsch@caltech.edu>, 10/15/2011.
# Purpose: given a big baroque output of exons from WormMart, get a simple nonredundant MERGED set of coordinates for exons from protein-coding genes in the genome.

use strict;
use warnings;

my $data_ref;

my $seq     = q{};
my $nt1     = 0;
my $nt2     = 0;
my $max     = 0;
my $flagged = 0;

while (my $input = <>) {
    chomp $input;

# Input fields:
#    Gene WB ID      Gene Public Name        Sequence Name (Transcript)      Transcript Type 
#                Chr Name        Strand  Start (bp)      End (bp)        
#                Exon Rank       Exon Start (bp) Exon End (bp)   
#                5' Intergenic Start (bp)        5' Intergenic End (bp)  3' Intergenic Start (bp)        3' Intergenic End (bp)

    if ( $input =~ /\A WBGene\d+ \t 
                       (?: [^\t]* \t){2}
                       coding \t 
                       (\S+) \t
                       (?: [^\t]* \t){4}
                       (\d+) \t (\d+) \t /xms ) { 
        $seq = $1;
        $nt1 = $2;
        $nt2 = $3;
        if ( $nt1 > $nt2 ) { 
            ($nt1, $nt2) = ($nt2, $nt1);
        }
        if ( $nt2 >= $max ) { 
            $max = $nt2;
        }
        # Note that I am assuming here that Perl will cope gracefully with 1-nt sites, for which $nt1 equals $nt2!
        foreach my $i ($nt1..$nt2) { 
            $data_ref->{'seq'}->{$seq}->{'flagged'}->{$i} = 1;
        }
    }
}

# Reset:
$nt1 = 0;
$nt2 = 0;

my @seqs = sort keys %{ $data_ref->{'seq'} };

foreach my $seq1 (@seqs) { 
    # Make the private loop variable public, so I can print it later if need be.
    $seq = $seq1;
    $flagged = 0;
    foreach my $j (1..$max) { 
        if (! $flagged) { 
            if ( exists $data_ref->{'seq'}->{$seq}->{'flagged'}->{$j} ) { 
                $nt1     = $j;
                $flagged = 1;
            }
            # Or, must be (! exists $data_ref->{'seq'}->{$seq}->{'flagged'}->{$j} ), 
            #     with no need to start a flagged zone, so do nothing.
        }
        # Otherwise, must be ($flagged):
        else {
            if ( exists $data_ref->{'seq'}->{$seq}->{'flagged'}->{$j} ) { 
                $nt2 = $j;
            }
            # Or, must be (! exists $data_ref->{'seq'}->{$seq}->{'flagged'}->{$j} ), 
            #     which means that flagged zone has just ended, so:
            else { 
                # Print out latest stored values for $seq, $nt1, $nt2:
                print "$seq\t$nt1\t$nt2\n";

                # Rezero everything:
                $nt1 = 0;
                $nt2 = 0;
                $flagged = 0;
            }
        }
    }
}

# Deal with "uncleared data when running off end of file" problem:
if ($flagged) {
    print "$seq\t$nt1\t$nt2\n";
}

