#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    # lcl|CM021144.1_cds_KAF1768189.1_1 [protein=hypothetical protein] [frame=2] [protein_id=KAF1768189.1] [location=join(<1..200,1085..1095)] [gbkey=CDS]
    if ( $input =~ /\A > /xms ) {
        if ( $input =~ /\A > .+ \[ protein_id = (\S+) \. \d+ \] .+ \z/xms ) {
            my $seqname = $1;
            my $header  = $input;
            $header  =~ s/\A[>]//;
            $header  = '>' . "$seqname  $header";
            $input   = $header;
        }
        else {
            die "Cannot parse FASTA header line: $input\n";
        }
    }
    print "$input\n";
}

