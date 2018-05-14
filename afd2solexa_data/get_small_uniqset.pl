#!/usr/bin/perl

# get_small_uniqset.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/18/2008.
# Purpose: with a GFF3 for *one* test gene, get just uniqs overlapping it (for parse-tests).

use strict;
use warnings;

my $input      = q{};
my $chr_test   = q{};
my $nt1_test   = q{};
my $nt2_test   = q{};
my $chromosome = 'I|II|III|IV|V|X';

while ($input = <>) { 
    chomp $input;

    # Typical tab-delimited input:
    # X  Coding_transcript  gene  7774799  7776622 [etc. etc.]

    if ( $input =~ /\A ( $chromosome ) 
                    \t Coding_transcript \t gene 
                    \t (\d+) 
                    \t (\d+) /xms ) { 
        ($chr_test, $nt1_test, $nt2_test) = ($1, $2, $3);
    }

    # Typical space-delimited input:
    # chrII 13068693 13068717 read 0 + - - 0,0,255

    if ( ($input) and ( $input =~ /\A chr($chromosome) \s+ (\d+) \s+ (\d+) /xms ) ) { 
        my $chr;
        my $nt1;
        my $nt2;
        ($chr, $nt1, $nt2) = ($1, $2, $3);
        if ( $chr eq $chr_test ) { 
            if ( fits_nts($nt1) or fits_nts($nt2) ) { 
                print "$input\n";
            }
        }
    }
}

sub fits_nts { 
    my $nt_to_check = $_[0];
    if (      ( $nt_to_check >= $nt1_test )
          and ( $nt_to_check <= $nt2_test ) ) { 
        return 1;
    }
}

