#!/usr/bin/env perl

# rename_seqs_for_phylip.pl -- Erich Schwarz <ems394@cornell.edu>, 11/23/2013.
# Purpose: convert names in a FASTA alignment into something much simpler that can go into a PHYLIP alignment; alternatively, make a synonym table, so that one can eventually take the PHYLIP-generated tree and put back full-length names.

use strict;
use warnings;

my $i = 0;

my $prefix = q{};
my $infile = q{};
my $output = q{};

if ( (! $ARGV[0]) or (! $ARGV[1]) or (! $ARGV[2]) ) {
    die "Format: rename_seqs_for_phylip.pl [prefix] [infile] [output type -- either 'fasta' or 'synonyms']\n";
}

$prefix = $ARGV[0];
$infile = $ARGV[1];
$output = $ARGV[2];

if ( $prefix !~ /\A \S+ \z/xms ) { 
    die "Can't parse prefix: $prefix\n";
}

if ( ( $output ne 'fasta' ) and ( $output ne 'synonyms' ) ) {
    die "Output must be either 'fasta' or 'synonyms'\n"
}

open my $INFILE, '<', $infile or die "Can't open infile $infile: $!";
while (my $input = <$INFILE>) { 
    chomp $input;
    if ( $input =~ /\A > (\S+) /xms ) { 
        my $orig_name = $1;
        $i++;
        my $synonym = $prefix . $i;
        print ">$synonym\n"            if ( $output eq 'fasta' );
        print "$synonym\t$orig_name\n" if ( $output eq 'synonyms' );
    }
    else {
        print "$input\n" if ( $output eq 'fasta' );
    }
}
close $INFILE or die "Can't close filehandle to infile $infile: $!";

