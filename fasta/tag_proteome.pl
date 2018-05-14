#!/usr/bin/perl

# tag_FASTA_names.pl
# Erich Schwarz, <emsch@its.caltech.edu>, 3/20/08
# Purpose: append tags like '_pcap' onto header names of a FASTA file.

use strict;
use warnings;

unless ( $#ARGV == 1 ) { 
    die 'Format: ./tag_FASTA_names.pl',
        '  [tag to append, e.g., \'_pcap\']',
        '  [input FASTA file]',
        "\n",
        ; 
}

my $tag = shift @ARGV;
if ($tag !~ / \A \w+ \z /xms) { 
     die "Cannot use tag \"$tag\".\n";
}
 
while (my $input = <>) { 
    chomp $input;
    if ($input !~ / \A > \S /xms) {
        print "$input\n";
    }
    elsif ($input =~ / \A > (\S+) (.*?) \z /xms) {
        my $input1 = $1;
        my $input2 = $2;
        print ">";
        print $input1;
        print $tag;
        print $input2;
        print "\n";
    }
    else { 
        die "Can't parse input line: $input\n";
    }
}

