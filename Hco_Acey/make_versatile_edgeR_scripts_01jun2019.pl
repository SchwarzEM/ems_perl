#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Getopt::Long;

use List::MoreUtils qw(uniq);

my $working_dir = q{};
my $input_data  = q{};
my @genotypes   = ();
my @comparisons = ();
my $suffix      = q{};

my $help;

GetOptions ( 'working_dir=s'    => \$working_dir,
             'input_data=s'     => \$input_data,
             'genotypes=s{,}'   => \@genotypes,
             'comparisons=s{,}' => \@comparisons,
             'suffix=s'         => \$suffix,
             'help'             => \$help,   );

if ( $help or (! $working_dir) or (! $input_data ) or (! @genotypes ) ) {
    die "Format: make_versatile_edgeR_scripts_01jun2019.pl \n",
        "    --working_dir|-w  [working directory for R]\n",
        "    --input_data|-i   [input data file for edgeR]\n",
        "    --genotypes|-g    [list of genotypes for input data replicate columns]\n",
        "    --comparisons|-c  [list of genotype pairs to be compared]\n",
        "    --suffix|-s       [suffix for output files]\n",
        "    --help|-h         [print this message]\n",
        ;
}

# We need a *non*-redundant genotype list, from which we can make comparison pairs and an input genotype text vector.

my @uniq_genotypes   = uniq(@genotypes);
my @quoted_genotypes = map { q{"} . $_ . q{"} } @genotypes;
my $genotype_vector  = join ",", @quoted_genotypes;

my @lines = ();

my $i = 0;

while (@comparisons) {
    my $geno1 = pop @comparisons;
    my $geno2 = pop @comparisons;

    my $line1 = '#' ." $geno1 vs. $geno2:";
    push @lines, $line1;

    # comp_XXX_vs_YYY <- exactTest(geno_comp, pair=c("XXX","YYY"))
    my $line2 = 'comp_' . $geno1 . '_vs_' . $geno2 . ' <- exactTest(geno_comp, pair=c' . "\(\"$geno1\",\"$geno2\"\)\)";
    push @lines, $line2;

    # mockdeg_XXX_vs_YYY <- topTags(comp_XXX_vs_YYY, n=Inf, p.value=1)
    my $line3 = 'mockdeg_' 
                 . $geno1 
                 . '_vs_' 
                 . $geno2 
                 . ' <- topTags(comp_' . $geno1 . '_vs_' . $geno2  . ', n=Inf, p.value=1)'
                 ;
    push @lines, $line3;

    # write.csv(as.data.frame(mockdeg_XXX_vs_YYY), file="XXX.vs.YYY_edgeR_exactTest_all.data_$suffix.csv")
    my $line4 = 'write.csv(as.data.frame(mockdeg_' 
                . $geno1 
                . '_vs_' 
                . $geno2 
                . '), file=' 
                . "\"$geno1.vs.$geno2" 
                . '_edgeR_exactTest_all.data_' 
                . $suffix 
                . ".csv\")"
                ;
    push @lines, $line4;
}

# Print a bunch of one-time stuff:

print '# Load edgeR:', "\n";
print 'library(edgeR)', "\n";
print '# Import data:', "\n";
print 'setwd("', $working_dir, '")', "\n";

# one very long line:
print 'edgeR_input_data';
print ' <- read.delim("', $input_data, '",row.names="Gene")', "\n";

# another very long line:
print 'genotypes <- factor(c(', $genotype_vector, '))', "\n";
print 'geno_comp <- DGEList(counts=edgeR_input_data,group=genotypes)', "\n";
print '# Check initial results:', "\n";
print 'geno_comp', "\n";
print '# Filter out weakly expressed genes using the following commands:', "\n";
print 'keep <- rowSums(cpm(geno_comp)>1) >= 3', "\n";
print 'geno_comp <- geno_comp[keep, , keep.lib.sizes=FALSE]', "\n";
print 'geno_comp$samples', "\n";
print '# Set up statistical data set:', "\n";
print 'geno_comp <- calcNormFactors(geno_comp)', "\n";
print 'geno_design <- model.matrix(~genotypes)', "\n";
print 'geno_comp <- estimateDisp(geno_comp,geno_design)', "\n";

# now, print out many stanzas of repetitive comparisons of different genotypes:

foreach my $repetitive_command (@lines) {
    print "$repetitive_command\n";
}

# some more one-time printing
print '# Sum up:', "\n";
print 'sessionInfo()', "\n";
print '# end it:', "\n";
print 'q()', "\n";

