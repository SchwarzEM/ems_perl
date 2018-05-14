#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

if (!@ARGV) {
    die "Usage: merge_columns.pl first_data_column.txt second_data_column.txt added_merge_text > [tab-linked two-column text]\n";
}

my $filenameA  = q{};
my $filenameB  = q{};
my $merge_text = q{};

$filenameA  = $ARGV[0] if $ARGV[0];
$filenameB  = $ARGV[1] if $ARGV[1];
$merge_text = $ARGV[2] if $ARGV[2];

if (! -r $filenameA ) {
    die "Could not read $filenameA\n";
}
if (! $filenameB) {
    die "Could not read $filenameB\n";
}

open my $FILEA, '<', $filenameA;
open my $FILEB, '<', $filenameB;

my $lineA = <$FILEA>;
my $lineB = <$FILEB>;

while(defined $lineA) {
        chomp $lineA;
	print "$lineA";
	$lineA = <$FILEA>;

        chomp $lineB;
	print "$merge_text\t$lineB\n";
	$lineB = <$FILEB>;
}
