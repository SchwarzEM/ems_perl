#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my @genotypes = ("WT_Mar2017", "vos2_Mar2017", "VOS2-GR_Untreat.4hr", "VOS2-GR_Treat.4hr");

my @lines = ();

my %have_compared = ();

my $i = 0;

foreach my $i (0..2) {
    foreach my $j ( ($i+1) .. 3 ) {
        my $geno1 = $genotypes[$i];
        my $geno2 = $genotypes[$j];

        my $line1 = '#' ." $geno1 vs. $geno2:";
        push @lines, $line1;

        # comp_XXX_vs_YYY <- exactTest(vos2_geno_comp, pair=c("XXX","YYY"))
        my $line2 = 'comp_' . $geno1 . '_vs_' . $geno2 . ' <- exactTest(vos2_geno_comp, pair=c' . "\(\"$geno1\",\"$geno2\"\)\)";
        push @lines, $line2;

        # mockdeg_XXX_vs_YYY <- topTags(comp_XXX_vs_YYY, n=Inf, p.value=1)
        my $line3 = 'mockdeg_' . $geno1 . '_vs_' . $geno2 . ' <- topTags(comp_' . $geno1 . '_vs_' . $geno2  . ', n=Inf, p.value=1)';
        push @lines, $line3;

        # write.csv(as.data.frame(mockdeg_XXX_vs_YYY), file="XXX.vs.YYY_edgeR_exactTest_all.data_2017.04.14.01.csv")
        my $line4 = 'write.csv(as.data.frame(mockdeg_' . $geno1 . '_vs_' . 
            $geno2 . '), file=' . "\"$geno1.vs.$geno2" . "_edgeR_exactTest_all.data_2017.04.14.01.csv\")";
        push @lines, $line4;
    }
}

# Print a bunch of one-time stuff:

print '# Load edgeR:', "\n";
print 'library(edgeR)', "\n";
print '# Import data:', "\n";
print 'setwd("/home/bioinformatics/VOS2_sequence_data")', "\n";

# one very long line:
print 'vos2_input_2017.04.14.01';
print ' <- read.delim("/home/bioinformatics/VOS2_sequence_data/vos2_edgeR.data_29mar2016.tsv.txt",row.names="Gene")', "\n";

# another very long line:
print 'genotypes <- factor(c("WT_Mar2017","WT_Mar2017","vos2_Mar2017","vos2_Mar2017",';
print '"VOS2-GR_Untreat.4hr","VOS2-GR_Untreat.4hr","VOS2-GR_Untreat.4hr","VOS2-GR_Treat.4hr","VOS2-GR_Treat.4hr","VOS2-GR_Treat.4hr"))', "\n";

print 'vos2_geno_comp <- DGEList(counts=vos2_input_2017.04.14.01,group=genotypes)', "\n";
print '# Check initial results:', "\n";
print 'vos2_geno_comp', "\n";
print '# Filter out weakly expressed genes using the following commands:', "\n";
print 'keep <- rowSums(cpm(vos2_geno_comp)>1) >= 3', "\n";
print 'vos2_geno_comp <- vos2_geno_comp[keep, , keep.lib.sizes=FALSE]', "\n";
print 'vos2_geno_comp$samples', "\n";
print '# Set up statistical data set:', "\n";
print 'vos2_geno_comp <- calcNormFactors(vos2_geno_comp)', "\n";
print 'vos2_geno_design <- model.matrix(~genotypes)', "\n";
print 'vos2_geno_comp <- estimateDisp(vos2_geno_comp,vos2_geno_design)', "\n";

# now, print out many stanzas of repetitive comparisons of different genotypes:

foreach my $repetitive_command (@lines) {
    print "$repetitive_command\n";
}

# some more one-time printing
print '# Sum up:', "\n";
print 'sessionInfo()', "\n";
print '# end it:', "\n";
print 'q()', "\n";

