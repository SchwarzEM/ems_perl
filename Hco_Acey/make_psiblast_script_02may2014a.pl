#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my @queries = qw(
    indiv_xenologs/Acey_s0010.g910.t3_40-156.fa
    indiv_xenologs/Acey_s0010.g910.t3_192-306.fa
    indiv_xenologs/Acey_s0010.g910.t3_335-449.fa
    indiv_xenologs/Acey_s0010.g910.t3_504-618.fa
    indiv_xenologs/Acey_s0010.g910.t3_678-792.fa

    indiv_xenologs/Acey_s0230.g2988.t2_36-112.fa
    indiv_xenologs/Acey_s0230.g2988.t2_145-259.fa
    indiv_xenologs/Acey_s0230.g2988.t2_311-424.fa
    indiv_xenologs/Acey_s0230.g2988.t2_482-595.fa
    indiv_xenologs/Acey_s0230.g2988.t2_652-766.fa

    indiv_xenologs/Acey_s0004.g1962.t1_57-145.fa
    indiv_xenologs/Acey_s0004.g2039.t1_3-109.fa
    indiv_xenologs/Acey_s0065.g3635.t1_51-168.fa
    indiv_xenologs/Acey_s0065.g3635.t1_197-274.fa

    indiv_xenologs/Acey_s0010.g996.t1_38-149.fa
    indiv_xenologs/Acey_s0212.g2239.t1_98-216.fa
    indiv_xenologs/Acey_s0517.g2812.t1.fa
);

my $header = '#!/bin/bash' . "\n\n";

foreach my $query (@queries) { 
    print $header if $header;
    $header = q{};
    my $stem = $query;
    $stem    =~ s/\.fa\z//;
    print "    psiblast -db dbs/metazoan_prots_14apr2014_orig.c_lectin_domains",
          " -query $query -num_threads 8 -evalue 1e-20 -inclusion_ethresh 1e-20 -num_iterations 20 -outfmt 6",
          " -seg yes -out $stem.psiblast_20x_1e-20.metazoan_prots_14apr2014_orig.c_lectin_domains.tsv.txt ;\n",
          ;
}
print "\n";

