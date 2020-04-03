#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A[>]/xms ) {
        # Sample header line:
        # >lcl|CM021144.1_cds_KAF1768189.1_1 [protein=hypothetical protein] [frame=2] [protein_id=KAF1768189.1] [location=join(<1..200,1085..1095)] [gbkey=CDS]
        if ( $input =~ /\A [>] (lcl\|\S+\.1_cds_ (\S+) \.1_\d+) (.*) \z/xms   ) {
             my $orig_id  = $1;
             my $short_id = $2;
             my $comments = $3;
             $input = '>' . "$short_id  $orig_id$comments"
        }
        else {
           warn "Cannot revise name in header: $input\n";
        }
    }
    print "$input\n";
}

