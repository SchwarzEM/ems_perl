#!/usr/bin/env perl

use strict ;
use warnings ;
use autodie ;

use List::MoreUtils qw(uniq);

my @genes = ();

while (my $input = <> ) {
    chomp $input;
    if ( $input =~ /\A [>] (\S+)\.t\d+ /xms ) {
        my $gene = $1;
        push @genes, $gene;
    }
    elsif ( $input =~ /\A [>] /xms ) {
        die "Cannot parse FASTA header line: $input\n";
    }
}

@genes = uniq(@genes);

if (@genes) {
    print "Gene\n";
    foreach my $gene (@genes) {
        print "$gene\n";
    }
}

