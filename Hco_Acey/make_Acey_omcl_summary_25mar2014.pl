#!/usr/bin/env perl

use strict;
use warnings;

my @key_spp = qw( ancylostoma 
                  elegans );

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
print "    rm Acey_omcl_summary_25mar2014.txt ;\n";
print "    touch Acey_omcl_summary_25mar2014.txt ;\n";
print "\n";

foreach my $key_species (@key_spp) { 
    my @other_spp = grep { $_ ne $key_species } @all_spp;
    foreach my $other_species (@other_spp) { 
        print "    echo \"Strict pairwise orthologs of $key_species vs. $other_species (other spp. completely variable):\" >> Acey_omcl_summary_25mar2014.txt ;\n";
        print '    /mnt/home/emsch/perl.svn/trunk/orthomcl_and_phylo/strict_vs_exp_omcls.pl',
              ' -i /mnt/home/emsch/work/Acey/ng_revision_mar2014/post_orthomcl/Acey_orthomcl_13spp_16mar2014.genes.txt',
              " -r $key_species $other_species",
              ' -a | /mnt/home/emsch/perl.svn/trunk/orthomcl_and_phylo/genes_from_taxon_in_omcl.pl -s -i -',
              " -t $key_species | grep 'Genes from taxon:' >> Acey_omcl_summary_25mar2014.txt ;\n",
              ;
        print "    echo >> Acey_omcl_summary_25mar2014.txt ;\n";
        print "\n";
    }
}

