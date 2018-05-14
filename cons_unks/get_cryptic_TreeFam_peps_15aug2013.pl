#!/usr/bin/env perl

use strict;
use warnings;

print '#!/bin/bash', "\n\n";

while (my $input = <>) { 
    chomp $input;
    print "    wget http://www.treefam.org/family/", $input, "/alignment ;\n";
    print "    mv alignment $input", "_raw.pep.fa ;\n";
    print "    uniform_fasta.pl -c '-' -i $input", "_raw.pep.fa | tag_FASTA_names.pl -i - -p $input", "_ > $input.pep.fa ;\n";
    print "    rm $input", "_raw.pep.fa ;\n";
    print "\n";
}

