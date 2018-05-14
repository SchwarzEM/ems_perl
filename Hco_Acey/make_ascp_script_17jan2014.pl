#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

print '#!/bin/bash', "\n\n";

while (my $input = <>) { 
    chomp $input;
    my $filename = basename($input);
    print "    /woldlab/rattus/lvol0/mus/home/schwarz/.aspera/connect/bin/ascp",
          " -i /sternlab/redivivus/data02/schwarz/H_contortus_genome/orig_reads_for_SRA/sra-8.ppk",
          " -QT -l600m -k1 $input",
          ' asp-sra@upload.ncbi.nlm.nih.gov:incoming ;',
          "\n",
          ;
    print "    e_ping -p done_upload.$filename.reads_to_SRA_17jan2014 ;\n";
    print "\n";
}

