#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use File::Basename;

print '#!/bin/bash', "\n\n";
print "    rm quick_report_2015.07.28.01.txt ;\n\n";

while (my $input = <>) {
    chomp $input;
    my $filename = basename($input);
    chomp $input;
    print "    echo \"Non-protein-coding expression data for ", 
          $filename, 
          ":\" >> quick_report_2015.07.28.01.txt ;\n",
          ;
    print "    cat $input | cut -f 1,6,10 | grep -v AT >> quick_report_2015.07.28.01.txt ;\n";
    print "    echo >> quick_report_2015.07.28.01.txt ;\n";
    print "\n";
}

