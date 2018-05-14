#!/usr/bin/env perl

# sieve_rscouts_w_blastx.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/24/2008.
# Purpose: filter RepeatScout hits if non-RT/transposase wormpep190 hit.

use strict;
use warnings;

my $blastx = $ARGV[0];
my $rep_fa = $ARGV[1];

my $print  = 1;
my %reject = ();

open my $BLASTX, '<', $blastx 
    or die "Can't open RScout-BlastX table $blastx: $!";

LOOP:
while ( my $input = <$BLASTX> ) { 
    chomp $input;
    if ( ( $input =~ / (R|r)everse \s transcriptase /xms)
          or ( $input =~ / (T|t)ransposase /xms) ) { 
        next LOOP;
    }
    if ( $input =~ / \A (\S+) .+ WBGene\d+ .+ \z /xms ) { 
        my $repeat = $1;
        $reject{$repeat} = 1;
    }
}
close $BLASTX 
    or die "Can't close filehandle for",
           " RScout-BlastX table $blastx: $!",
           ;

open my $REP_FA, '<', $rep_fa
    or die "Can't open FASTA $rep_fa for filtering: $!";
while ( my $input = <$REP_FA> ) { 
    chomp $input;
    if ( $input =~ / \A > (\S+) /xms ) { 
        my $repeat = $1;
        if ( $reject{$repeat} ) { 
            $print = 0;
        }
        else { 
            $print = 1;
            print "$input\n";
        }
    }
    else { 
        if ($print) { 
            print "$input\n";
        }
    }
}
close $REP_FA or die "Can't close filehandle to FASTA $rep_fa: $!";

