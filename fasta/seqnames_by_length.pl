#!/usr/bin/env perl

# seqnames_by_length.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/22/2008.
# Purpose: from FASTA, list seq. names and sizes in descending size order.

use strict;
use warnings;

my $name    = q{};
my $length  = 0;
my %seq2len = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > (\S+) /xms ) {
        my $new_name = $1;
        if ($name) { 
            $seq2len{$name} = $length;
            $length = 0;
        }
        $name = $new_name;
        $new_name = q{};
    }
    if ( ( $input !~ /\A >/xms ) and ( $input =~ /[a-zA-Z]/xms ) ) { 
        # Formerly was: $input =~ s/\A\s+|\s+\z//; -- what was I thinking?
        $input =~ s/\s//g;
        if ( $input =~ /[^a-zA-Z]/xms ) {
            die "Illegal input line: $input\n";
        }
        $length += length($input);
    }
}

if ($name) {
    $seq2len{$name} = $length;
    $length = 0;
}

my @sorted_keys = sort { $seq2len{$b} <=> $seq2len{$a} } keys %seq2len;
foreach my $seq ( @sorted_keys ) { 
    print "$seq\t$seq2len{$seq}\n";
}

