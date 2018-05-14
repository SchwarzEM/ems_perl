#!/usr/bin/env perl

use strict;
use warnings;

my %num2spp = ( '9606' => 'human',
                '10090' => 'mouse',
                '7955' => 'zebrafish',
                '7227' => 'melanogaster',
                '6239' => 'elegans',
                '44689' => 'dictyostelium',
                '559292' => 'cerevisiae',
                '284812' => 'pombe',
                '3702' => 'arabidopsis',
);

my @id_numbers = sort keys %num2spp;

print '#!/bin/bash -login', "\n";
print '#PBS -l walltime=006:00:00', "\n";
print '#PBS -l nodes=1:ppn=1', "\n";
print '#PBS -l mem=10gb', "\n";
print '#PBS -N job_wrangle_pfam_10aug2013', "\n";
print '#PBS -q main', "\n";
print '#PBS -M ems394@cornell.edu', "\n";
print '#PBS -m abe', "\n";
print '#PBS -A ged-intel11', "\n";
print '#PBS -r n', "\n";
print '#PBS -V', "\n";
print "cd /mnt/home/emsch/work/DUF_R01/pfam/precomps ;\n";
print "gunzip *.gz ;\n";

foreach my $id (@id_numbers) { 
    my $new_id = $num2spp{$id} . '_' . $id;
    print "mv $id.tsv $new_id.tsv ;\n";
}

