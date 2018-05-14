#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $family = $ARGV[0];
my $input  = $ARGV[1];
my $print  = 0;

if ( (! $family) or (! $input) ) { 
    die "Format: extract_one_treefam_family_09mar2014.pl [one TreeFam family accession] [input monster TreeFam FASTA with accessions] > [raw FASTA for the one family]\n";
}

open my $INPUT, '<', $input;
while (my $input = <$INPUT>) { 
    chomp $input;
    if ( $input =~ /\A (TF\d+) \s* \z/xms ) { 
        my $new_family = $1;
        if ( $new_family eq $family ) { 
            $print = 1;
        }
        else { 
            $print = 0;
        }
    }
    elsif ($print) { 
        print "$input\n" if ( $input !~ /\A [#]END /xms );
    }
}
close $INPUT;

