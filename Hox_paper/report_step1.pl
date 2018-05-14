#!/usr/bin/perl

# report_step1.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/10/2007.
# Purpose: step 1 of churning out Tables 1-2 for 2007 Hox paper.
# Uses old inputs to avoid breaking stuff downstream; but has simple tab-delimited output.

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

    # Number of residues (nt, aa) in $seq:
    my $seqlen = q{};

    # Print extra line before starting new fosmid's data:
    if ($seq =~ /\. tfa \z/xms ) {
        print "\n";
    }

    print "$seq\t";

    # Weird regex-magic from Perl Cookbook 2.16, p. 84:
    $seqlen = commify($fa_lens{$seq});

    # Append residue type based solely on suffix -- kludge!
    if ($seq =~ /\. tfa \z/xms ) {
        $seqlen .= " nt";
    }
    if ($seq !~ /\. tfa \z/xms ) {
        $seqlen .=  " aa";
    }

    print "$seqlen\n";
}

# Perl Cookbook 2.16, p. 84 -- baffling, but it works:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

