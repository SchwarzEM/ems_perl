#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

my $output = 'indiv_reads_fastq_count_28jun2012a.txt';

print '#!/bin/bash', "\n\n";
print "   echo >  $output;\n\n";

while (my $input = <>) { 
    chomp $input;
    my $basename = basename $input;
    print "   echo \"fastq count of $basename\" >> $output;\n";
    print "   count_simple_fastq_residues.pl $input | grep Total >> $output;\n";
    print "   echo >> $output;\n";
    print "\n";
}
    
print "program_done_e-ping.pl -p counting_raw_indiv_LC_fastqs ;\n\n";


