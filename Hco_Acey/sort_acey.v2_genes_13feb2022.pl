#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my $header = "Gene\n";

while (my $input = <>) {
    chomp $input;
    # Sample input name: Acey_s0001.v2.g17940
    if ( $input =~ /\A Acey_s(\d+)\.v2\.g(\d+) \z/xms ) { 
        my $scaf_no  = $1;
        my $locus_no = $2;
        $scaf_no     =~ s/\A[0]+//g;
        $scaf_no     = (1000000 * $scaf_no);
        my $gene_no  = ($scaf_no + $locus_no);
        $data_ref->{'index'}->{$gene_no}->{'gene_id'} = $input;
    }
    elsif ( $input !~ /\A Gene \z/xms ) {
        die "Cannot process input line: $input\n";
    }
}

my @indexes = sort { $a <=> $b } keys %{ $data_ref->{'index'} };
foreach my $index (@indexes) {
    print $header if $header;
    $header = q{};
    my $gene_id = $data_ref->{'index'}->{$index}->{'gene_id'};
    print "$gene_id\n";
}


