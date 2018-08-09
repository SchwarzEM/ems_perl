#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A [>] (.*) \z/xms ) {
        my $header_text = $1;
        if ( $header_text =~ /\A (\S+) \b (.*) \z/xms ) { 
            my $orig_name = $1;
            my $comments  = $2;
            if ( $comments =~ / locus=(\S+)/xms ) {
                my $new_name = $1;
                print ">$new_name  $orig_name$comments\n";
            }
            else {
                $orig_name =~ s/[a-z]\.\d+\z//;
                $orig_name =~ s/[a-z]\z//;
                print ">$orig_name$comments\n";
            }
       }
    }
    elsif ( $input =~ /\A [>] /xms ) {
        die "Cannot parse input: $input\n";
    }
    else {
        print "$input\n";
    }
}

