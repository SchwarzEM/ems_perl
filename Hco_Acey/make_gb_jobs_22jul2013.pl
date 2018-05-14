#!/usr/bin/env perl

use strict;
use warnings;

print '#!/bin/bash', "\n\n";
print '    PATH="$PATH:/home/schwarz/perl.svn/trunk/ncbi" ;', "\n\n";

foreach my $i (1..7) {
    my $input_file = 'Hco_v4_coding.v7.part.' . $i . '.fa';

    print q{    faNs2ncbi.pl -m 1 -c -i }, 
          $input_file, 
          q{ | perl -ne '$_ =~ s/Hco_v4_scaf_/Hco_v4_contig_/g; print $_;' > }, 
          "Hco_v4_coding.v7.part.$i.contigs_19jul2013.fsa ;\n",
          ;

    print q{    faNs2agp.pl -m 1 -i }, 
          $input_file, 
          q{ | perl -ne 'my $input = $_; chomp $input; if ( $input =~ /\A ( Hco_v4_scaf .+? ) Hco_v4_scaf ( .+ ) \z /xms )},
          " {",
          q{ print "$1"; print "Hco_v4_contig"; print "$2\n"; },
          "}",
          q{ else },
          "{",
          q{ print "$input\n" },
          "}' > Hco_v4_coding.v7.part.$i.agp ;\n",
          ;

    print q{    tbl2asn -t Hco_17jul2013.sbt -i },
          $input_file,
          q{  -a s -j "[organism=Haemonchus contortus] [strain=McMaster] [tech=wgs]" -V v -M n},
          " -Z discrep_Hco.v7.part.$i.22jul2013.txt ;\n",
          ;
    print "  gzip -9 Hco_v4_coding.v7.part.$i.contigs_19jul2013.fsa Hco_v4_coding.v7.part.$i.agp Hco_v4_coding.v7.part.$i.sqn ;\n";
    print "\n";
}
    
print "    e_ping -p done_gb_jobs_22jul2013;\n\n";
