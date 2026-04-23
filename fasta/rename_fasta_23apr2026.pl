#!/usr/bin/env perl

# Purpose: given a renaming table and a FASTA file, rename its sequences while preserving its original FASTA text as headers.

use strict;
use warnings;
use autodie;

my $fasta = q{};
my $table = q{};

$fasta = $ARGV[0] if  $ARGV[0];
$table = $ARGV[1] if  $ARGV[1];

my %rename = ();

if ( (! $fasta ) or (! $table ) ) {
    die "Format: rename_fasta_23apr2026.pl [FASTA] [renaming table] > [renamed FASTA]\n";
}

open my $TABLE, '<', $table;
while ( my $input = <$TABLE> ) {
    chomp $input;
    my $orig = q{};
    my $rev  = q{};
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        $orig = $1;
        $rev  = $2;
        if ( exists $rename{$orig} ) {
            die "From renaming table $table, redundant renaming of $orig: $rename{$orig} and $rev\n";
        }
        $rename{$orig} = $rev;
    }
    else {
        die "From renaming table $table, cannot parse: $input\n";
    }
}
close $TABLE;

open my $FASTA, '<', $fasta;
while ( my $input = <$FASTA> ) {
    chomp $input;

    my $orig = q{};
    my $rev  = q{};
    my $head = q{};

    if ( ( $input !~ /\A [>] /xms ) and ( $input =~ /\S/xms ) ) {
        print "$input\n";
    }
    elsif ( $input =~ /\A [>] ((\S+) .*) \z/xms ) {
        $head = $1;
        $orig = $2;
        if (! exists $rename{$orig} ) {
            warn "In FASTA $fasta, could not rename $orig in: $input\n";
            $rev = $orig;
        }
        else {
            $rev = $rename{$orig};
        }
        print '>' . "$rev  $head\n";
    }
}
close $FASTA;
