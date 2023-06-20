#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A \S+ \t \S+ \t gene \t+ .+ \t [^\t]* ID=([^;\s]+)/xms ) {
        my $gene_id = $1;
        $input = "$input;note=$gene_id";
    }
    elsif ( $input =~ /\A \S+ \t \S+ \t mRNA \t+ .+ \t [^\t]* ID=([^;\s]+)/xms ) {
        my $tx_id = $1;
        $input = "$input;note=$tx_id";
    }
    elsif ( $input =~ /\A .+ \t [^\t]* ID=[^;\s]+ .* transcript_id=([^;\s]+) /xms ) {
        my $prot_id = $1;
        $prot_id =~ tr/[a-z]/[A-Z]/;
        $input = "$input;note=$prot_id";
    }
    print "$input\n";
}
