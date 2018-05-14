#!/usr/bin/env perl

# cds2wbfullgenes.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/19/2008.
# Purpose: given an incoming CDS list stream and a defined wormpep, emit full ID stream.

use strict;
use warnings;

my %cds2fullname = ();
my %possible_cds = ();
my %seen         = ();

# Typical wormpep190 lines:
# >4R79.1b        CE39659 WBGene00003525  locus:nas-6 ..
# >AC7.3  CE07653 WBGene00014997

while (my $input = <>) { 
    chomp $input;
    if ($input =~ /\A > (\S+) \s+ CE\d+ \s+ (WBGene\d+) \s+ (\S+) \s+  /xms) { 
        my $cds         = $1;
        my $gene        = $2;
        my $maybe_locus = $3;
        my $locus       = q{};
        $cds =~ s/[a-z]\z//;
        my @id_tags = ($gene, $cds);
        if ( $maybe_locus =~ /\A locus : (\S+) \z/xms ) { 
            $locus = $1;
            push @id_tags, $locus;
        }
        my $full_id = join '|', @id_tags;
        $cds2fullname{$cds} = $full_id;
    }
    elsif ( $input =~ /\A (\S+) (\s.*) \z /xms ) { 
        my $cds          = $1;
        my $rest_of_line = $2;
        if ( $cds2fullname{$cds} ) { 
            print $cds2fullname{$cds},
                  $rest_of_line,
                  "\n",
                  ;
        }
    }
}

