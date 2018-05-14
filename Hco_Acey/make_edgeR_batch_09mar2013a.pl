#!/usr/bin/env perl

use strict;
use warnings;

my @pairwises = ( q{1,2},
                  q{1,3},
                  q{2,3},
                  q{3,4},
                  q{4,5},
                  q{5,6},
                  q{6,15},
                  q{8,7},
                  q{10,9},
                  q{12,11},
                  q{14,13}, );

print q{library(edgeR)}, "\n";
print q{setwd("/sternlab/redivivus/data02/schwarz/Acey_genomics/post_meltdown/edgeR_2013.03.09.01")}, "\n";
print q{Acey_hkeep_counts         <- read.delim("Acey_pme_expected_count.rounded.5plus_reads.01mar2013.406genes.txt",row.names="Gene")}, "\n";
print q{Acey_hkeep_counts.DGEList <- DGEList(counts=Acey_hkeep_counts)}, "\n";
print q{Acey_hkeep_counts.DGEList <- calcNormFactors(Acey_hkeep_counts.DGEList)}, "\n";
print q{Acey_hkeep_counts.DGEList <- estimateCommonDisp(Acey_hkeep_counts.DGEList)}, "\n";
print q{Acey_counts         <- read.delim("Acey_pme_expected_count.rounded.5plus_reads.01mar2013.txt",row.names="Gene")}, "\n";
print q{Acey_grouping       <- factor(c(1:15))}, "\n";
print q{Acey_counts.DGEList <- DGEList(counts=Acey_counts,group=Acey_grouping)}, "\n";
print q{Acey_counts.DGEList <- calcNormFactors(Acey_counts.DGEList)}, "\n";
print q{Acey_counts.DGEList$common.dispersion <- Acey_hkeep_counts.DGEList$common.dispersion}, "\n";

foreach my $pairing (@pairwises) { 
    my $pair_label = $pairing;
    $pair_label    =~ s/,/v/;
    print 'Acey_counts.', $pair_label, '.DGEList.exactTest <- exactTest(Acey_counts.DGEList, pair=c(', $pairing, '))', "\n", ;
    print 'Acey_counts.', $pair_label, '.DGEList.exactTest.topTags <- topTags(Acey_counts.', $pair_label, '.DGEList.exactTest, n=Inf)', "\n", ;
    print 'write.table(Acey_counts.', $pair_label, '.DGEList.exactTest.topTags$table, file="Acey_counts.', $pair_label, '.DGEList.exactTest.topTags.txt")', "\n", ;
}

print "q()\n";


