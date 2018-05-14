#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use File::Basename;

my @infiles = ();

while (my $input = <>) {
    chomp $input;
    push @infiles, $input;
}
        
print '#!/bin/bash -login', "\n";
print '#PBS -l walltime=024:00:00', "\n";
print '#PBS -l nodes=1:ppn=1', "\n";
print '#PBS -l mem=32gb', "\n";
print "#PBS -N job_compress_and_md5sum_sra_inputs_18feb2016.sh\n";
print '#PBS -q main', "\n";
print '#PBS -M ems394@cornell.edu', "\n";
print '#PBS -m abe', "\n";
print '#PBS -A ged', "\n";
print '#PBS -r n', "\n";
print '#PBS -V', "\n";

foreach my $infile (@infiles) {
    my $dir = dirname($infile);
    my $file = basename($infile);
    print "cd $dir ;\n";
    print "gzip $file ;\n";
    print "md5sum $file.gz > $file.gz.md5sum ;\n";
}

