#!/usr/bin/env perl

use strict;
use warnings;

print '#!/bin/bash', "\n\n";
print "    mkdir sqn_products_23jul2013a ;\n\n";

foreach my $i (1..7) {
    my $input_file = 'Hco_v4_coding.v7.part.' . $i . '.contigs_19jul2013.fsa';
    print q{    tbl2asn -t Hco_17jul2013.sbt -i },
          $input_file,
          q{  -a s -j "[organism=Haemonchus contortus] [strain=McMaster] [tech=wgs]" -V v -M n},
          " -Z discrep_Hco.v7.part.$i.23jul2013a.txt ;\n",
          ;
    print "\n";
}

print "    mv *.sqn sqn_products_23jul2013a ;\n";
print "    gzip -9 sqn_products_23jul2013a/*.sqn ;\n\n";
print "    e_ping -p done_gb_jobs_23jul2013a;\n\n";

