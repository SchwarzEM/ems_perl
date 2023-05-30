#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A ((\S+)_fimo)\/fimo.gff \z/xms ) {
        my $dir  = $1;
        my $stem = $2;
        my $outfile = $dir . q{/} . $stem . '_fimo.bed'; 
        print "cat $input | gff2bed > $outfile ;\n";
    }
    else {
        die "Cannot parse input: $input\n";
    }
}

