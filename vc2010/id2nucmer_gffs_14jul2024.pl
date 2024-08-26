#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <> ) {
    chomp $input;
    if ( $input !~ /\A[#]/xms ) {
        if ( $input =~ /\A (\S+ \t \S+) \t \S+ \t (.*) \z/xms ) {
            my $text1 = $1;
            my $text2 = $2;
            $input = "$text1\tregion\t$text2";
        }
        else {
            die "Cannot parse third column in: $input\n";
        }

        if ( $input =~ /\A (\S+) \t .*? \t \S+ \t (\d+) \t (\d+) \t /xms ) {
            my $seq = $1;
            my $nt1 = $2;
            my $nt2 = $3;
            my $ident = "ID=$seq:$nt1-$nt2";
            $input =~ s/\tTarget=\S+[ ]\d+[ ]\d+[;]/\t/;
            $input =~ s/\t([^\t]*)\z/\t$ident;$1/;
            $input =~ s/Name=/Previous=/;
        }
        else {
            die "Cannot parse seq. and coords. in: $input\n";
        }
    }
    print "$input\n";
}

