#!/usr/bin/perl

# 100res_lengths.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/19/2007.
# Purpose: feeble hack for listing FASTA-formatted sequences only if >=100 residues; *assumes* not junk.

use strict;
use warnings;

my $fasta   = q{};
my %fa_lens = ();

while (my $input = <>) { 
    chomp $input;
    if ($input =~ / \A > (\S+) /xms) { 
        $fasta = $1;
        if ($fa_lens{$fasta}) { 
            die "Multiple entries for $fasta\n";
        }
        $fa_lens{$fasta} = 0;
    }
    elsif ( ( $input =~ /\S/ ) and ($fasta) ) { 
        $input =~ s/\s//g;
        $fa_lens{$fasta} += length($input);
    }
}

foreach my $seq (sort keys %fa_lens) { 
    if ($fa_lens{$seq} >= 100) { 
        print "$seq",
              "\t",
              "$fa_lens{$seq}",
              "\n",
              ;
    }
}

