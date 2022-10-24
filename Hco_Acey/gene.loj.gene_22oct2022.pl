#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $data_ref;

my $header = "Gene\tMapped_genes";
my @genes1 = ();

while (my $input = <>) {
    if ($input =~ /ID[=](\S+?)[;] .+ ID[=](\S+)/xms ) {
        my $gene1 = $1;
        my $gene2 = $2;
        push @genes1, $gene1;
        $data_ref->{'gene1'}->{$gene1}->{gene2}->{$gene2} = 1;
    }
}

@genes1 = uniq @genes1;
foreach my $gene1 (@genes1) {
     my @genes2 = sort keys %{ $data_ref->{'gene1'}->{$gene1}->{gene2} };
     my $gene2_text = join '; ', @genes2;
     if ($header) {
         print "$header\n";
         $header = 0;
     }
     print "$gene1\t$gene2_text\n";
}
