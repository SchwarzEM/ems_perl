#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <> ) {
    chomp $input;
    my $output = q{};
    if ( $input !~ /\A[#]/xms ) {
        # Sample input line:
        # chrI    94910   147548
        # Intended output line:
        # chrI    94910   147548  .       .       .       maf-convert     region  .       ID=chrI:94911-147548;prev_id_no=5

        if ( $input =~ /\A (\S+) \t (\d+) \t (\d+) \z/xms ) {
            my $seq   = $1;
            my $nt1   = $2;
            my $nt2   = $3;

            # Because this is a BED rather than a GFF, there is a weird -1 effect on the starting nt which must be reversed:
            $nt1++;

            my $ident = "ID=$seq:$nt1-$nt2";            
            $output = "$input\t.\t.\t.\tmaf-convert\tregion\t.\t$ident";       
        }
        else {
            die "Cannot parse/revise: $input\n";
        }
    }
    print "$output\n";
}

