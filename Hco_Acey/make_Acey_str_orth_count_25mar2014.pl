#!/usr/bin/env perl

use strict;
use warnings;

my @key_spp = qw( ancylostoma 
                  elegans );

# Note that this list explicitly is *not* including ascaris_jex -- we are using the later genome exclusively for this analysis.
my @all_spp = qw( ancylostoma
                  ascaris
                  ascaris_jex
                  briggsae 
                  brugia
                  b_xylophilus
                  dirofilaria
                  elegans
                  haemonchus
                  m_hapla
                  necator 
                  pristionchus
                  trichinella );

print '#!/bin/bash', "\n\n";

foreach my $key_species (@key_spp) { 
    print "    rm $key_species.genes.w_strict_orths.orig.txt ;\n";
    print "    touch $key_species.genes.w_strict_orths.orig.txt ;\n"; 
    my @other_spp = grep { $_ ne $key_species } @all_spp;
    foreach my $other_species (@other_spp) { 
        print '    /mnt/home/emsch/perl.svn/trunk/orthomcl_and_phylo/strict_vs_exp_omcls.pl',
              ' -i /mnt/home/emsch/work/Acey/ng_revision_mar2014/post_orthomcl/Acey_orthomcl_13spp_16mar2014.genes.txt',
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

