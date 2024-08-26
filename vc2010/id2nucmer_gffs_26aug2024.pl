#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <> ) {
    chomp $input;
    if ( $input !~ /\A[#]/xms ) {
        # Sample input line:
        # chrI    maf-convert     region  3533    37088   .       .       .       ID=1
        if ( $input =~ /\A (\S+) \t maf-convert \t region \t (\d+) \t (\d+) \t .* \t ID=(\d+) \z/xms ) {
            my $seq   = $1;
            my $nt1   = $2;
            my $nt2   = $3;
            my $id_no = $4;

            my $ident = "ID=$seq:$nt1-$nt2";            
            $input =~ s/\t[^\t]*\z/\t$ident;prev_id_no=$id_no/;
        }
        else {
            die "Cannot parse/revise: $input\n";
        }
    }
    print "$input\n";
}

