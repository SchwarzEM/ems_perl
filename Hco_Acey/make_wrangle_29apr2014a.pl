#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = '#!/bin/bash' . "\n\n";

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (nonmet_amiDs_blastp_\S+).fa.txt \z/xms ) { 
        my $stem   = $1;
        my $outlist = $stem . '.top50list.txt';
        my $output  = $stem . '.top50list.fa';
        print $header if $header;
        $header = q{};
        print "    grep '>' $input | perl -ne '",
              q{ $input = $_; if ( $input =~ /\A > (\S+) /xms ) { print "$1\n"; } '},
              " | head --lines=50 > $outlist ;\n",
              ;
        print "    extract_fasta_subset.pl -l $outlist -f $input | niceify_ncbi_proteomes_29apr2014.pl > $output ;\n";
        print "\n";
    }
}

