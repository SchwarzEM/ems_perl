#!/usr/bin/env perl

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A \S+ \/ ([^\s\/]+) \.genes\.results \z /xms ) {
        my $tag = $1;
        my $short_tag = $tag;
        # /Hco_v4_cDNA.Adult.female.lib_11445 to Adult.female 
        $short_tag =~ s/\AHco_v4_cDNA\.//;
        $short_tag =~ s/\.lib_\d+//;
        $short_tag =~ s/L2-both/L2/;
        $short_tag =~ s/L3_lib_Melb/L3/;
        print "    cut -f 1,9 $input",
              ' | perl -ne \' s/\bpme_TPM\b/',
              $short_tag,
              '_TPM/; print \' | perl -ne \' s/\bgene_id\b/Gene/; print \'',
              " > $tag.genes.pme_TPMs ;\n",
              ;
    }
}

