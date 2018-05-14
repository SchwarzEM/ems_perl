#!/usr/bin/env perl

# add_fasta_suffixes.pl -- Erich Schwarz <ems394@cornell.edu>, 12/4/2015.
# Purpose: given a FASTA file with redundant sequence names, add suffixes to make the names non-redundant.

use strict;
use warnings;
use autodie;

my $infile = $ARGV[0];

my %count  = ();
my %suffix = ();

# count names
open my $INFILE, '<', $infile;
while (my $input = <$INFILE>) {
    if ( $input =~ /\A > (\S+) /xms ) { 
        my $name = $1;
        $count{$name}++;
    }
}
close $INFILE;

# given a full namecount, give them suffixes *if* they have redundant names.
open $INFILE, '<', $infile;
while (my $input = <$INFILE>) {
    chomp $input;
    if ( $input =~ /\A > (\S+) /xms ) {
        my $name = $1;
        if ( $count{$name} >= 2 ) {
            $suffix{$name}++;
            my $revname = $name . '_mod' . $suffix{$name};
            $input =~ s/\A[>]$name/>$revname/;
        }
    }
    print "$input\n";
}
close $INFILE;

