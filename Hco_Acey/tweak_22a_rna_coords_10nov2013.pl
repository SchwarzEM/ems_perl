#!/usr/bin/env perl

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

my $coords      = $ARGV[0];
my $start_tweak = $ARGV[1];
my $stop_tweak  = $ARGV[2];

if (    (! looks_like_number($start_tweak) ) or ( int($start_tweak) != $start_tweak ) 
     or (! looks_like_number($stop_tweak) )  or ( int($stop_tweak)  != $stop_tweak ) ) { 
    die "Start and stop tweaks both need to be integers.\n";
}

open my $COORDS, '<', $coords or die "Can't open DNA FASTA file with embedded header coordinates, $coords: $!";
while (my $input = <$COORDS>) { 
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        if ( $input =~ /\A > ([^\|\s]+) \| (\d+) \| (\d+) \s* \z/xms ) {
            my $seq      = $1;
            my $start_nt = $2;
            my $end_nt   = $3;

            $start_nt = $start_nt + $start_tweak;
            $end_nt   = $end_nt   + $stop_tweak;

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

