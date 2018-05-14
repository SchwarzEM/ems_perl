#!/usr/bin/perl

# uniform_fasta.pl -- Erich Schwarz <emsch@its.caltech.edu>, 6/30/2006.
# Purpose: filter a FASTA file (with >=1 seqs.) into simple alph.-ord. headers and ALL-CAPS 60-char/line.

use strict;
use warnings;

my $input_line   = "";
my $output_line  = "";
my @output_lines = ();
my $seq_name     = "";
my %sequences    = ();

if ($#ARGV > 0) { die "This script is merging multiple files into one output!\n"; } 

while (<>) { 
    chomp ($input_line = $_);
    if ($input_line =~ /^>(\S+)/) { 
        $seq_name = $1; 
        $sequences{$seq_name} = "";
    }
    elsif ($input_line =~ /[a-zA-Z]/) { 
        $input_line =~ s/[^a-zA-Z]//g;
        $input_line =~ tr/[a-z]/[A-Z]/;
        $sequences{$seq_name} .= $input_line;
    }
}

foreach $seq_name (sort keys %sequences) { 
    print ">$seq_name\n";
    @output_lines = unpack("a60" x (length($sequences{$seq_name})/60 + 1), $sequences{$seq_name});
    foreach $output_line (@output_lines) { 
          print "$output_line\n";
    }
}

