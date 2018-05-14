#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

my $output = 'indiv_mapped_genes_count_28jun2012c.txt';

print '#!/bin/bash', "\n\n";
print "    echo > $output ;\n\n";

while (my $input = <>) { 
    chomp $input;
    my $basename = basename $input;
    print "    echo \"Count of mapped genes in $basename:\" >> $output ;\n";
    print "    cat $input ", q{| perl -ne 'm/\A(\S+)\s+\S+\s+(\S+)/; $wbgene = $1; $rpkm = $2; if ( $rpkm > 0 ) { print "$wbgene\n"; } ' | sort | uniq | wc -l }, ">> $output ;\n", ; 
    print "    echo >> $output ;\n";
    print "\n";
}
    
print "    program_done_e-ping.pl -p counting_mapped_indiv_LC_genes ;\n\n";

