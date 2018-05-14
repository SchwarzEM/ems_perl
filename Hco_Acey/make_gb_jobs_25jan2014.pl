#!/usr/bin/env perl

use strict;
use warnings;

print '#!/bin/bash', "\n\n";
print '    PATH="$PATH:/home/schwarz/perl.svn/trunk/ncbi" ;', "\n\n";

foreach my $i (1..4) {
    my $input_file = 'Acey_2013.11.30.genDNA.ncbi_scaf.rev1.part.' . $i . '.fa';

    print q{    faNs2ncbi.pl -m 1 -c -i }, 
          $input_file, 
          q{ | perl -ne '$_ =~ s/_scaf/_contig/g; print $_;' > }, 
          "Acey_2013.11.30.genDNA.ncbi_scaf.rev1.part.$i.contigs_25jan2014.fsa ;\n",
          ;

    print q{    tbl2asn -t Acey_15jan2014.sbt -i },
          "Acey_2013.11.30.genDNA.ncbi_scaf.rev1.part.$i.contigs_25jan2014.fsa",
          q{ -a s -j "[organism=Ancylostoma ceylanicum] [tech=wgs]" -V v -M n},
          " -Z discrep_Acey_2013.11.30.genDNA.ncbi_scaf.rev1.part.$i.contigs_25jan2014.txt ;\n",
          ;
    print "  gzip -9 Acey_2013.11.30.genDNA.ncbi_scaf.rev1.part.$i.contigs_25jan2014.fsa",
          " Acey_2013.11.30.genDNA.ncbi_scaf.rev1.part.$i.contigs_25jan2014.sqn ;\n",
          ;
    print "\n";
}
    
print "    e_ping -p done_gb_jobs_25jan2014;\n\n";
