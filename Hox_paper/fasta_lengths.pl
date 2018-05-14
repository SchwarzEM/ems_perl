#!/usr/bin/perl

# fasta_lengths.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/19/2007.
# Purpose: quick hack to get lengths of FASTA-formatted sequences; *assumes* they're not junk.

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
    print "$seq",
          "\t",
          "$fa_lens{$seq}",
          "\n",
          ;
}

