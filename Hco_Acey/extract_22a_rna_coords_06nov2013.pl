#!/usr/bin/env perl

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

my $flanks = $ARGV[0];
my $coords = $ARGV[1];

if ( (! looks_like_number($flanks) ) or ( int($flanks) != $flanks ) or ( $flanks < 1 ) ) { 
    die "Flank size $flanks needs to be a positive integer.\n";
}

open my $COORDS, '<', $coords or die "Can't open DNA FASTA file with embedded header coordinates, $coords: $!";
while (my $input = <$COORDS>) { 
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        if ( $input =~ /\A > ([^\|\s]+) \| (\d+) \| (\d+) \s* \z/xms ) {
            my $seq      = $1;
            my $start_nt = $2;
            my $end_nt   = $3;

            $start_nt = $start_nt + $flanks;
            $end_nt   = $end_nt   - $flanks;

            if ( $start_nt >= $end_nt ) {
                die "In DNA FASTA file with embedded header coordinates, $coords, can't make sense of coordinates: $input\n";
            }
            print "$seq\t$start_nt\t$end_nt\n";
        }
        else { 
            die "In DNA FASTA file with embedded header coordinates, $coords, can't parse header line: $input\n";
        }
    }
}
close $COORDS or die "Can't close fileheader to DNA FASTA file with embedded header coordinates, $coords: $!";

