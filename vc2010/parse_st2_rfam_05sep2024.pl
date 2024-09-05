#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

while (my $input = <> ) {
    chomp $input;
    # Note that we are trimming ".1", ".2" etc. from StringTie2 transcript names here to get 'gene' names.
    if ( $input =~ /\A ([^\t]+) \t ([^\t]+) \t ([^\t]+) \.\d+ \t ([^\t]+) \z/xms ) {
        my $mot_name  = $1;
        my $mot_acc   = $2;
        my $query_seq = $3;
        my $mot_desc  = $4;
        my $motif = "$mot_acc|$mot_name|\"$mot_desc\"";
        $data_ref->{'query_seq'}->{$query_seq}->{'motif'}->{$motif} = 1;
    }
    else {
        die "Cannot parse input: $input\n";
    }
}

my @query_seqs = sort keys %{ $data_ref->{'query_seq'} };
foreach my $query_seq (@query_seqs) {
    my @motifs = sort keys %{ $data_ref->{'query_seq'}->{$query_seq}->{'motif'} };
    foreach my $motif (@motifs) {
         print "$query_seq\t$motif\n";
    }
}

