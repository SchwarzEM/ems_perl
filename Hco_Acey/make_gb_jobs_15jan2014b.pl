#!/usr/bin/env perl

use strict;
use warnings;

print '#!/bin/bash', "\n\n";

foreach my $i (1..4) {
    my $input_file = 'Acey_2013.11.30.genDNA.ncbi_scaf.part.' . $i . '.contigs_15jan2014.fsa';
    print q{    tbl2asn -t Acey_15jan2014.sbt -i },
          $input_file,
          q{  -a s -j "[organism=Ancylostoma ceylanicum] [tech=wgs]" -V v -M n},
          " -Z discrep_Acey_2013.11.30.genDNA.ncbi_scaf.part.$i.contigs_15jan2014.txt ;\n",
          ;
    print "\n";
}

print "    gzip -9 *.sqn ;\n\n";
print '    gzip -9 Acey_2013.11.30.genDNA.ncbi_scaf.part.*', ".contigs_15jan2014.fsa ;\n\n";
print "    e_ping -p done_gb_jobs_15jan2014b;\n\n";

