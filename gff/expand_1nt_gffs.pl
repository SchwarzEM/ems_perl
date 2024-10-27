#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while ( my $input = <> ) {
    chomp $input;
    if ( $input =~ /\A [# ]/xms ) {
        print "$input\n";
    }
    # Sample input:
    # Sherm_2022.07.24.06     maf-convert     region  9550    10446   .       .       .       ID=1
    elsif ( $input =~ /\A ( (\S+) \t [^\t]* \t [^\t]*) \t (\d+) \t (\d+) \t (.*) \z /xms ) {
        my $text1 = $1;
        my $seq   = $2;
        my $a     = $3;
        my $b     = $4;
        my $text2 = $5;

        if ( $a > $b ) {
            die "$a is larger than $b in: $input\n";
        }

        foreach my $i ($a..$b) {
            my $id = 'ID=' . "$seq.nt-$i";
            # change ID=x to ID=seq.nt
            if ( $text2 =~ /ID=\S+/xms ) {
                $text2 =~ s/ID=\S+/$id/g;
            }
            elsif ( $text2 =~ /\A (.+ \t) ([^\t*]) \x/xms ) {
                my $subtext2a = $1;
                my $subtext2b = $2;
                $text2 = $subtext2a . "$id;" . $subtext2b;
            }
            print "$text1\t$i\t$i\t$text2\n";
        }
    }
    else {
        die "Cannot parse input: $input\n";
    }
}
