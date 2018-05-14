#!/usr/bin/env perl

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A \.\. \/ Hco_(\S+)_rsem_runs \/ ([^\s\/]+) \.genes\.results \z /xms ) {
        my $version = $1;
        my $tag     = $2;
        my $short_tag = $tag;

        # /Hco_v4_cDNA.Adult.female.lib_11445 to Adult.female for short tag; Hco_v4_cDNA to Hco_v4_Pasi or Hco_v4_aug for long tag.

        $tag       =~ s/\AHco_v4_cDNA\./Hco_v4_$version./;

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

