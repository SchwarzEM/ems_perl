#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A[>](\S+)/xms ) {
        my $seqid = $1;
        $input = '>' . "$seqid [organism=Steinernema hermaphroditum] [strain=PS9179]";
        if ( $seqid =~ /\A \S+ _ [A-Z]+\z/xms ) {
            $input = "$input [chromosome=$seqid]";
        }
    }
    print "$input\n";
}
