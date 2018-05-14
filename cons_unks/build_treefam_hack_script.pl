#!/usr/bin/env perl

use strict;
use warnings;
        
print '#!/bin/bash -login', "\n";
print '#PBS -l walltime=003:00:00', "\n";
print '#PBS -l nodes=1:ppn=1', "\n";
print '#PBS -l mem=10gb', "\n";
print '#PBS -N ', "job_build_treefam_hack_14aug2013\n";
print '#PBS -q main', "\n";
print '#PBS -M ems394@cornell.edu', "\n";
print '#PBS -m abe', "\n";
print '#PBS -A ged-intel11', "\n";
print '#PBS -r n', "\n";
print '#PBS -V', "\n";
print "cd /mnt/home/emsch/work/DUF_R01/treefam/hacked_TreeFams ;\n";

while (my $input = <>) { 
    chomp $input;
    print 'wget http://www.treefam.org/family/', $input, "/tree/newick ;\n";
    print "mv newick $input.newick.txt ;\n";
}

