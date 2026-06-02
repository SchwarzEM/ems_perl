#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

while (my $input = <> ) {
    chomp $input;
    if ( $input =~ /\A Gene \t/xms ) {
        print "Gene\tSex-expr-bias\n";
    }
    elsif ( $input =~ /\A (\S+) \t ([^\t]*) \t ([^\t]*) \z/xms ) {
        my $gene = $1;
        my $fc   = $2;
        my $fdr  = $3;
        if ( ( looks_like_number($fc) ) and ( looks_like_number($fdr) ) and ( $fdr <= 0.001 ) ) {
            if ( $fc >= 1 ) {
                print "$gene\tMale_bias\n";
            }
            elsif ( $fc <= -1 ) {
                print "$gene\tFem_bias\n";
            }
            else {
                print "$gene\tNo_bias\n";
            }
        }
        else {
            print "$gene\tNo_bias\n";
        }
    }
    else {
        die "Cannot parse: $input\n";
    }
}
