#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $orig_table = $ARGV[0];
my $pval_table = $ARGV[1];

my %comp2pval  = ();
my $comparison = q{};

# Load this value in to make header easier later:
$comp2pval{'Comparison'} = 'p_value';

open my $PVAL_TABLE, '<', $pval_table;
while (my $input = <$PVAL_TABLE>) {
    chomp $input;
    # Sample input line:
    # > # Comparison: "atml1-3.vs.lgo-2.up";
    if ( $input =~ /\A [>] \s+ [#] \s+ Comparison[:] \s+ \" (\S+) \" [;] \s+ /xms ) {
        $comparison = $1;
    }
    # Sample input line:
    # number of successes = 23, number of trials = 89, p-value < 2.2e-16
    elsif ( $input =~ /\A number \s+ of \s+ successes \s+ [=] \s+ \d+, 
                          \s+ number \s+ of \s+ trials \s+ [=] \s+ \d+, 
                          \s+ p-value \s+ ([(?:<|=)] \s \S+) \s* \z/xms ) {
        my $p_value = $1;
        # Delete '=', but keep (and trim down) '<'.
        $p_value =~ s/\A[=]\s+//;
        $p_value =~ s/\A[<]\s+/</;
        if ( $comparison =~ /\S/xms ) {
            if ( exists $comp2pval{$comparison} ) { 
                die "Redundant p-values for $comparison: $comp2pval{$comparison} vs. $p_value\n";
            }
            if ( $p_value eq 'TRUE' ) { 
                $p_value = 1;
            }
            $comp2pval{$comparison} = $p_value;
            $comparison             = q{};
        }
    }
}
close $PVAL_TABLE;

open my $ORIG_TABLE, '<', $orig_table;
while (my $input = <$ORIG_TABLE>) {   
    chomp $input;
    if ( $input =~ /\A (\S+) (\t .+) \z/xms ) {
        my $comparison = $1;
        my $text       = $2;
        if (! exists $comp2pval{$comparison} ) {
            die "Cannot find p-value for comparison (\"$comparison\") in: $input\n";
        }
        print "$input\t$comp2pval{$comparison}\n";
    }
}
close $ORIG_TABLE;


