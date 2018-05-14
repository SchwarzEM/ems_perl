#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = '#!/bin/bash' . "\n\n";
my $tail   = q{};

while (my $input = <>) {
    chomp $input;
    my $tag = q{};
    if ( $input =~ /rsem_ (\S+) _2015/xms ) { 
        $tag = $1;
    }
    else {
        die "Can't parse input: $input\n";
    }
    if ($header) {
        print $header;
        $header = q{};
        $tail   = "\n";
    }   
    print "    cut -f 1,10 $input",
          ' | perl -ne \' $input = $_; $input =~ s/gene_id/Gene/; $input =~ s/pme_TPM/',
          $tag,
          '_pme_TPM/; print $input; \' > partial_tables/',
          $tag,
          "_pme_TPM.tsv.txt ;\n",
          ;
}

print $tail;

