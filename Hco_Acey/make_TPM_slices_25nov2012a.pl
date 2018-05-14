#!/usr/bin/env perl

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A \S+ \/ ([^\s\/]+) \.genes\.results \z /xms ) {
        my $tag = $1;
        print "    cut -f 1,9 $input",
              ' | perl -ne \' s/\bpme_TPM\b/',
              $tag,
              # From the following line, omit: q{ | perl -ne \' s/\b0.00\b/0.01/; print \' }
              '_TPM/; print \' | perl -ne \' s/\bgene_id\b/Gene/; print \'',
              " > $tag.genes.pme_TPMs ;\n",
              ;
    }
}

