#!/usr/bin/env perl

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (\S+\.fa) \.orig_revcomp \z/xms ) {
        my $new_name = $1;
        if (-e $new_name) { 
            die "Try revcomping $input to something that doesn't already exist, like $new_name\n";
        }
        print "    make_revcomp_seqs.pl -s -f $input > $new_name ;\n";
    }
    else { 
        die "Can't parse input: $input\n";
    }
}

