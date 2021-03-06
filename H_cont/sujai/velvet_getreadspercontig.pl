#!/usr/bin/env perl

# velvet_getreadspercontig.pl -- Sujai Kumar <sujai.kumar@ed.ac.uk> -- 7/10/2010.

use strict;
use warnings;

use Tie::File;

die
"Usage: velvet_getreadspercontig.pl {<velvet_directory>} {<file_with_contig_list>}\n"
  unless @ARGV == 2;

my $directory  = shift @ARGV;
my $contiglist = shift @ARGV;

open CONTIGS, "<$contiglist" or die "No file of contigs found\n";
my %contiglist;
while (<CONTIGS>) { chomp; $contiglist{$_}++; }
close CONTIGS;

open GRAPH, "<$directory/LastGraph" or die "$directory/LastGraph not found\n";

my %reads;

while (<GRAPH>) {
  NR: if (/^NR\t-?(\d+)/) {
        my $contig_id = $1;
        next unless exists $contiglist{$contig_id};
        while (<GRAPH>) {
            if   (/^(\d+)/) { $reads{$1}++; }
            else            { goto NR }
        }
    }
}

open SEQUENCES, "<$directory/Sequences"
  or die "$directory/Sequences not found\n";
while (<SEQUENCES>) {
    if ( /^>.*\s(\d+)\s+\d+$/ and exists $reads{$1} ) {
        print $_;
        my $seq = <SEQUENCES>;
        print $seq;
    }
}
