#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <> ) {
    chomp $input;
    my $seqname = q{};
    my $output = $input;
    if ( ( $output !~ /\A[#]/xms ) and ( $output =~ /\A (\S+) /xms ) ) {
        $seqname = $1;
        $output =~ s/=g/=$seqname.g/g; 
    }
    print "$output\n";
}

