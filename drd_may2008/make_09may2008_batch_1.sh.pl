#!/usr/bin/env perl

# make_02may2008_batch_1.sh.pl -- Erich Schwarz, 5/2/08.
# Purpose: make a simple-minded but serviceable batch shell script for DRD psi-blast.

# Read in directory contents:

use strict;
use warnings;

my $long_script_name = '/home/schwarz/taygeta'
                       . '/science/perl_scripts_svn'
                       . '/trunk/integr8_uniprot'
                       . '/extract_uniprot_subset.pl'
                       ;

# Hard-coded -- kludgey, but should work.

# Only change from 5/5: 1e-09 hits to 1e-12 hits.
my $select_list_name   = '../DRD.psiblast.all_integr8.30x.1e-12.hits.txt';
my $select_list_output = $select_list_name . '_integr8.dat';

my %dat_file_names = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ / \A \s* ( \S+\.dat) \.gz \z/xms ) { 
        my $dat_name = $1;
        $dat_file_names{$dat_name} = 1;
    }
}

print '#!/bin/bash',                "\n",
                                    "\n",
      'cd integr8_proteomes_gz ;',  "\n",
      "rm $select_list_output ;",   "\n",
                                    "\n",
      ;

foreach my $dat_file_name (sort keys %dat_file_names) { 
    print "zcat $dat_file_name.gz > $dat_file_name ;\n";
    print "$long_script_name  $select_list_name  $dat_file_name >> $select_list_output ;\n";
    print "rm $dat_file_name ;\n\n";
}

