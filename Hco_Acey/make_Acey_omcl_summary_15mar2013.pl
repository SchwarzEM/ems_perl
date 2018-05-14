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
print "    rm Acey_omcl_summary_15mar2013.txt ;\n";
print "    touch Acey_omcl_summary_15mar2013.txt ;\n";
print "\n";

foreach my $key_species (@key_spp) { 
    my @other_spp = grep { $_ ne $key_species } @all_spp;
    foreach my $other_species (@other_spp) { 
        print "    echo \"Strict pairwise orthologs of $key_species vs. $other_species (other spp. completely variable):\" >> Acey_omcl_summary_15mar2013.txt ;\n";
        print '    /home/schwarz/perl.svn/trunk/orthomcl/strict_vs_exp_omcls.pl',
              ' -i /sternlab/redivivus/data02/schwarz/Acey_genomics/augustus_preds_24oct2012/annots_1/Acey_orthomcl.10spp_27oct2012.genes.txt',
              " -r $key_species $other_species",
              ' -a | /home/schwarz/perl.svn/trunk/orthomcl/genes_from_taxon_in_omcl.pl -s -i -',
              " -t $key_species | grep 'Genes from taxon:' >> Acey_omcl_summary_15mar2013.txt ;\n",
              ;
        print "    echo >> Acey_omcl_summary_15mar2013.txt ;\n";
        print "\n";
    }
}

