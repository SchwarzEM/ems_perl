#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+ \t \S+ \t gene \t .+) \t (\S+) \s* \z/xms ) {
        my $lead_text = $1;
        my $id        = $2;
        $input = "$lead_text\tID=$id";
    }
    if ( $input =~ /\A (\S+ \t \S+ \t transcript \t .+) \t (\S+) \s* \z/xms ) {
        my $lead_text = $1;
        my $id        = $2;
        my $gene_id   = $id;
        $gene_id      =~ s/\.t\d+\z//;
        $input = "$lead_text\tID=$id; gene_id \"$gene_id\"";
    }
    print "$input\n";
}

