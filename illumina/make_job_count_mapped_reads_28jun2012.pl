#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

my $output = 'indiv_mapped_reads_count_28jun2012b.txt';

print '#!/bin/bash', "\n\n";
print "    echo > $output ;\n\n";

while (my $input = <>) { 
    chomp $input;
    my $basename = basename $input;
    print "    echo \"Count of mapped reads in $basename:\" >> $output ;\n";
    print "    cat $input ", q{| perl -ne 'm/\A(\S+)/; print "$1\n"; ' | sort | uniq | wc -l }, ">> $output ;\n", ; 
    print "    echo >> $output ;\n";
    print "\n";
}
    
print "    program_done_e-ping.pl -p counting_mapped_indiv_LC_fastqs ;\n\n";

