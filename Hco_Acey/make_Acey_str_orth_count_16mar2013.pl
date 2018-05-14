#!/usr/bin/env perl

use strict;
use warnings;

my @key_spp = qw( ancylostoma 
                  elegans );

my @all_spp = qw( ancylostoma
                  ascaris
                  b_xylophilus
                  briggsae
                  brugia
                  dirofilaria
                  elegans
                  m_hapla
                  pristionchus
                  trichinella );

print '#!/bin/bash', "\n\n";

foreach my $key_species (@key_spp) { 
    print "    rm $key_species.genes.w_strict_orths.orig.txt ;\n";
    print "    touch $key_species.genes.w_strict_orths.orig.txt ;\n"; 
    my @other_spp = grep { $_ ne $key_species } @all_spp;
    foreach my $other_species (@other_spp) { 
        print '    /home/schwarz/perl.svn/trunk/orthomcl/strict_vs_exp_omcls.pl',
              ' -i /sternlab/redivivus/data02/schwarz/Acey_genomics/augustus_preds_24oct2012/annots_1/Acey_orthomcl.10spp_27oct2012.genes.txt',
              " -r $key_species $other_species",
              ' -a | perl -ne \' while ( /(\S+\(',
              $key_species,
             '\))/g ) { print "$1\n"; } \' | sort | uniq',
             " >> $key_species.genes.w_strict_orths.orig.txt ;\n",
             ;
    }
    print "    sort $key_species.genes.w_strict_orths.orig.txt | uniq > $key_species.genes.w_strict_orths.txt ;\n";
    print "    wc -l $key_species.genes.w_strict_orths.txt > $key_species.genes.w_strict_orths.count.txt ;\n";
    print "\n";
}

