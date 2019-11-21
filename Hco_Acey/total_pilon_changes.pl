#!/usr/bin/env perl

# total_pilon_changes.pl -- Erich Schwarz <ems394@cornell.edu>, 11/22/2019.
# Purpose: given a change log from Pilon, get a count of all changes (i.e., both all deletions and all insertions).

use strict;
use warnings;
use autodie;

my $deletions  = 0;
my $insertions = 0;
my $total      = 0;
my $net        = 0;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A .+ \s+ (\S+) \s+ (\S+) /xms ) {
        my $prev_seq   = $1;
        my $rev_seq    = $2;
        my $prev_count = 0;
        my $rev_count  = 0;
        my $diff       = 0;

        $prev_seq =~ s/[^ACGT]//g;
        $rev_seq  =~ s/[^ACGT]//g;

        $prev_count = length($prev_seq);
        $rev_count  = length($rev_seq);
        $diff       = $rev_count - $prev_count;

        if ( $diff > 0 ) {
            $insertions = $insertions + $diff;
        }
        elsif ( $diff < 0 ) {
            $deletions = $deletions - $diff;
        }
    }
    else {
        die "Cannot parse input line: $input\n";
    }
}

$total = $insertions + $deletions;
$net   = $insertions - $deletions;

$insertions = commify($insertions);
$deletions  = commify($deletions);
$total      = commify($total);
$net        = commify($net);

print "Deletions:      $deletions nt\n";
print "Insertions:     $insertions nt\n";
print "Total changes:  $total nt\n";
print "Net changes:    $net nt\n";

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}


