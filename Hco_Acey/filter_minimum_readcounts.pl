#!/usr/bin/env perl

# filter_readcounts.pl -- Erich Schwarz <ems394@cornell.edu>, 6/1/2019.
# Purpose: given readcount table and defined minimum readcount, delete every gene/transcript row that fails to achieve that minimum.

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

my $threshold  = 0 ;
my $read_table = q{};

my $header = 1;
my $print  = 0;

$threshold  = $ARGV[0] if $ARGV[0];
$read_table = $ARGV[1] if $ARGV[1];

if ((! $threshold) and (! $read_table) ) {
    die "Format: filter_minimum_readcounts.pl [minimum readcount threshold] [readcount table] > [filtered readcount table]\n";
}
if (! looks_like_number($threshold) ) {
    die "First argument needs to be a numerical minimum readcount threshold.\n";
}
if ( $threshold < 0 ) {
    die "Minimum readcount threshold must be at least 0.\n";
}

open my $READS, '<', $read_table;
while (my $input = <$READS>) {
    chomp $input;
    # Print the very first (header) line of the table without question, one time.
    if ($header) {
        print "$input\n";
        $header = 0;
    }
    # Filter all subsequent data lines.
    else { 
        $print = 0;
        if ( $input =~ /\A \S+ \t (.+) \z/xms ) {
            my $read_line = $1;
            my @count_vals = split /\t/, $read_line;
            foreach my $count (@count_vals) {
                if (! looks_like_number($count) ) {
                    die "From reads table $read_table, cannot parse readcounts from: $input\n";
                }
                elsif ( $count >= $threshold ) {
                    $print = 1;
                }
            }
            if ($print) {
                print "$input\n";
            }
        }
        else {
           die "From reads table $read_table, cannot parse data line: $input\n";
        }
    }
}
close $READS;
