#!/usr/bin/env perl

# make_02may2008_batch_1.sh.pl -- Erich Schwarz, 5/2/08.
# Purpose: make a simple-minded but serviceable batch shell script for DRD psi-blast.

# Read in directory contents:

use strict;
use warnings;

my $longperlname = '/home/schwarz/taygeta'
                   . '/science/perl_scripts_svn'
                   . '/trunk/integr8_uniprot'
                   . '/uniprot_to_clean_tfa.pl';
my %dat_names = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ / \A \s* ( \S+\.dat) \.gz \z/xms ) { 
        my $dat_name = $1;
        $dat_names{$dat_name} = 1;
    }
}

print '#!/bin/bash',                                  "\n",
                                                      "\n",
      'cd integr8_proteomes_gz ;',                    "\n",
      'rm ../integr8_proteomes/raw_all_integr8.fa ;', "\n",
                                                      "\n",
      ;

foreach my $dat_name (sort keys %dat_names) { 
    print "zcat $dat_name.gz > $dat_name ;\n";
    print "$longperlname $dat_name >> ../integr8_proteomes/raw_all_integr8.fa ;\n";
    print "rm $dat_name ;\n\n";
}

