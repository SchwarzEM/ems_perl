#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $infile  = q{};
my $species = q{};
my $unique  = q{};

$infile  = $ARGV[0] if $ARGV[0];
$species = $ARGV[1] if $ARGV[1];
$unique  = $ARGV[2] if $ARGV[2];

if ( (! $infile ) or (! $species ) ) {
    die "Format: orthogroups2genes_15nov2023.pl [infile] [species] [optionally, 'unique'] > gene list (optionally only unique orthologs)\n";
}

open my $INFILE, '<', $infile;
while ( my $input = <$INFILE> ) {
    chomp $input;
    while ( $input =~ / (\S+) \($species\) /xmsg ) {
         my $ortholog = $1;
         if ( (! $unique ) or ( $input !~ / $species .+ $species /xms ) ) {
             print "$ortholog\n";
         }
    }
    if ( $input !~ / \S+ \($species\) /xms ) {
        die "Cannot parse input: $input\n";
    }
}
close $INFILE;

