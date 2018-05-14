#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A [>] ((FL81_\d+)\-R[A-Z]) \s* \z/xms ) { 
        my $prot = $1;
        my $gene = $2;
        $prot    = 'Crem_' . $prot;
        $gene    = 'Crem_' . $gene;
        print ">$prot    gene:$gene\n";
    }
    elsif ( $input =~ /\A [>] /xms ) {
        die "Cannot process header text: $input\n";
    }
    else {
        print "$input\n";
    }
}

