#!/usr/bin/env perl

use strict;
use warnings;

print '#!/bin/bash', "\n\n";
print '    PATH="$PATH:/home/schwarz/perl.svn/trunk/ncbi" ;', "\n\n";

foreach my $i (1..4) {
    my $input_file = 'Acey_2013.11.30.genDNA.ncbi_scaf.part.' . $i . '.fa';

    print q{    faNs2ncbi.pl -m 1 -c -i }, 
          $input_file, 
          q{ | perl -ne '$_ =~ s/_scaf/_contig/g; print $_;' > }, 
          "Acey_2013.11.30.genDNA.ncbi_scaf.part.$i.contigs_15jan2014.fsa ;\n",
          ;

    print q{    faNs2agp.pl -m 1 -i }, 
          $input_file, 
          q{ | perl -ne 'my $input = $_; chomp $input; if ( $input =~ /\A ( Acey_s\d+_scaf \s .+? ) (Acey_s\d+) _scaf ( .+ ) \z /xms )},
          " {",
          q{ print "$1"; print "$2"; print "_contig"; print "$3\n"; },
          "}",
          q{ else },
          "{",
          q{ print "$input\n" },
          "}' > Acey_2013.11.30.genDNA.ncbi_scaf.part.$i.agp ;\n",
          ;

    print q{    tbl2asn -t Acey_15jan2014.sbt -i },
          $input_file,
          q{  -a s -j "[organism=Ancylostoma ceylanicum] [tech=wgs]" -V v -M n},
          " -Z discrep_Acey_2013.11.30.genDNA.ncbi_scaf.part.$i.15jan2014.txt ;\n",
          ;
    print "  gzip -9 Acey_2013.11.30.genDNA.ncbi_scaf.part.$i.contigs_15jan2014.fsa",
          " Acey_2013.11.30.genDNA.ncbi_scaf.part.$i.agp",
          " Acey_2013.11.30.genDNA.ncbi_scaf.part.$i.sqn ;\n",
          ;
    print "\n";
}
    
print "    e_ping -p done_gb_jobs_15jan2014;\n\n";
