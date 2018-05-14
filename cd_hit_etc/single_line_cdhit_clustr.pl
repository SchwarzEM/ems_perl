#!/usr/bin/env perl

# single_line_cdhit_clustr.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/10/2010.
# Purpose: convert *.clstr outputs of CD-HIT into single-line summaries.

use strict;
use warnings;

my $read_line = 0;
my $seqname   = q{};
my $keyseq    = q{};
my %clustseqs = ();
my $synonyms  = q{};

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        $synonyms = q{};
        if (%clustseqs) { 
            $synonyms = join "\t", (sort keys %clustseqs);
            $synonyms = "\t" . $synonyms;
        }
        if ($read_line) { 
            print "$keyseq";
            print "$synonyms" if $synonyms;
            print "\n";
        }
        $seqname   = q{};
        $keyseq    = q{};
        %clustseqs = ();
        $synonyms = q{};
    }
    elsif ( $input =~ /\A [^>]* > (\S+) \.\.\. (.*) /xms ) {
        $seqname = $1;
        $read_line = 1;
        if ( $input =~ / \* \s* \z /xms ) { 
            $keyseq = $seqname;        
        }
        else { 
            $clustseqs{$seqname} = 1;
        }
    }
}

print "$keyseq";
print "$synonyms" if $synonyms;
print "\n";

