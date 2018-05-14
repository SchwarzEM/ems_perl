#!/usr/bin/env perl

use strict;
use warnings;

my $seqnames = $ARGV[0];
my $wormpep  = $ARGV[1];

my %bad_seqs = ();
my %bad_peps = ();

open my $SEQNAMES, '<', $seqnames or die "Can't open sequence names file $seqnames: $!";
while (my $input = <$SEQNAMES>) {
    chomp $input;
    $bad_seqs{$input} = 1;
}
close $SEQNAMES or die "Can't close filehandle to sequence names file $seqnames: $!";

open my $WORMPEP, '<', $wormpep or die "Can't open wormpep file $wormpep: $!";
while (my $input = <$WORMPEP>) {
    chomp $input;
    if ( $input =~ /\A > ([^\.\s]+ \. \d+ [a-z]) \s /xms ) { 
        my $pep = $1;
        my $seq = $pep;
        $seq    =~ s/[a-z]\z//;
        if ( $bad_seqs{$seq} ) {
            $bad_peps{$pep} = 1;
        }
    }
    elsif ( $input =~ /\A > ([^\.\s]+ \. \d+) \s /xms ) {
        my $pep = $1;
        if ( $bad_seqs{$pep} ) {
            $bad_peps{$pep} = 1;
        }
    }
    elsif ( $input =~ /\A > /xms ) { 
        die "Can't parse input header line: $input\n";
    }
}
close $WORMPEP or die "Can't close filehandle to wormpep file $wormpep: $!";

my @rejects = sort keys %bad_peps;

foreach my $badpep (@rejects) { 
    print "$badpep\n";
}
