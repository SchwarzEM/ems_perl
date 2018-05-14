#!/usr/bin/env perl

use strict;
use warnings;

print '#!/bin/bash', "\n\n";

foreach my $i (1..11) {
    # Give the $i two digits, for '01' to 11'.
    $i = sprintf ("%02u", $i);
    my $input_file = 'Acey_2013.11.30.ncbi.mod2.cDNA.part.' . $i . '.fsa';
    print q{    tbl2asn -t Acey_23jan2014.sbt -i },
          $input_file,
          q{  -a s -j "[organism=Ancylostoma ceylanicum] [moltype=transcribed_RNA] [tech=TSA]" -V v -M n},
          " -Z discrep_Acey_2013.11.30.ncbi.mod2.cDNA.part.$i.18jan2014.txt ;\n",
          ;
    print "    gzip -9 Acey_2013.11.30.ncbi.mod2.cDNA.part.$i.fsa",
          " Acey_2013.11.30.ncbi.mod2.cDNA.part.$i.sqn ;\n",
          ;
    print "\n";
}

print "    e_ping -p done_gb_jobs_23jan2014;\n\n";

