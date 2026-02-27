#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $infile = q{};
$infile    = $ARGV[0] if $ARGV[0];

if (! $infile ) {
    die "Format: top_bitscore_outfmt6_blast_26feb2026.pl [outfmt6 BLAST tabular output] > [only topmost bitscores per query]\n";
}

my $data_ref;

open my $INFILE, '<', $infile;
while ( my $data = <$INFILE> ) {
    chomp $data;
    if ( $data =~ /\A (\S+) \t .+ \t (\S+) \z/xms ) {
        my $seqid = $1;
        my $score = $2;
        if ( exists $data_ref->{'seqid'}->{$seqid}->{'max'} ) {
            my $max = $data_ref->{'seqid'}->{$seqid}->{'max'};
            if ( $score > $max ) {
                $data_ref->{'seqid'}->{$seqid}->{'max'} = $score;
                $data_ref->{'seqid'}->{$seqid}->{'data'} = $data;
            }
        }
        elsif (! exists $data_ref->{'seqid'}->{$seqid}->{'max'} ) {
            $data_ref->{'seqid'}->{$seqid}->{'max'} = $score;
            $data_ref->{'seqid'}->{$seqid}->{'data'} = $data;
        }
        else {
            die "Cannot assess maximum score of $seqid\n"
        }
    }
    else {
        die "From infile $infile, cannot parse: $data\n";
    }
}
close $INFILE;

my @seqids = sort keys %{ $data_ref->{'seqid'} };
foreach my $seqid2 (@seqids) {
    if ( exists $data_ref->{'seqid'}->{$seqid2}->{'data'} ) {
        my $data2 = $data_ref->{'seqid'}->{$seqid2}->{'data'};
        print "$data2\n";
    }
    else {
        die "Cannot map sequence $seqid2 to top-bitscore-selected data\n";
    }
}
