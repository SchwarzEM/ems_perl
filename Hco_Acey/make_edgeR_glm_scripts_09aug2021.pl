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
    die "Format: make_edgeR_glm_scripts_09aug2021.pl  [arguments] > [edgeR glm / pairwise-comparisons R script]\n",
        "    --working_dir|-w  [working directory for R]\n",
        "    --input_data|-i   [input data file for edgeR]\n",
        "    --genotypes|-g    [list of genotypes for input data replicate columns]\n",
        "    --comparisons|-c  [list of genotype pairs to be compared]\n",
        "    --suffix|-s       [suffix for output files]\n",
        "    --help|-h         [print this message]\n",
        ;
}

my @quoted_genotypes = map { q{"} . $_ . q{"} } @genotypes;
my $genotype_vector  = join ",", @quoted_genotypes;

my @lines = ();

my $i = 0;

while (@comparisons) {
    my $geno1 = shift @comparisons;
    my $geno2 = shift @comparisons;

    my $line1 = '#' ." $geno1 vs. $geno2:";
    push @lines, $line1;

    # comp_XXX_vs_YYY <- makeContrasts(XXX - YYY, levels=geno_design)
    # qlf_XXX_vs_YYY <- glmQLFTest(geno_comp_fit, contrast=comp_XXX_vs_YYY)
    # If XXX has stronger expression than YYY, it will get a positive logFC value.

    my $line2 = 'comp_' . $geno1 . '_vs_' . $geno2 . ' <- makeContrasts' . "\($geno2 - $geno1" . ', levels=geno_design)';
    push @lines, $line2;

    my $line3 = 'qlf_' . $geno1 . '_vs_' . $geno2 . ' <- glmQLFTest(geno_comp_fit, contrast=comp_' . $geno1 . '_vs_' . $geno2 . ')';
    push @lines, $line3;

    # deg_XXX_vs_YYY <- topTags(comp_XXX_vs_YYY, n=Inf, p.value=1)
    my $line4 = 'deg_' 
                 . $geno1 
                 . '_vs_' 
                 . $geno2 
                 . ' <- topTags(qlf_' . $geno1 . '_vs_' . $geno2  . ', n=Inf, p.value=1)'
                 ;
    push @lines, $line4;

    # write.table(as.data.frame(deg_XXX_vs_YYY), file="XXX.vs.YYY_edgeR_glm.pairwise_$suffix.csv", sep = "\t")
    my $line5 = 'write.table(as.data.frame(deg_' 
                . $geno1 
                . '_vs_' 
                . $geno2 
                . '), file=' 
                . "\"$geno1.vs.$geno2" 
                . '_edgeR_glm.pairwise_' 
                . $suffix 
                . ".orig.tsv.txt\""
                . ', sep="\t", col.names=NA)'
                ;
    push @lines, $line5;
}

# Print a bunch of one-time stuff:

print '# Load edgeR:', "\n";
print 'library(edgeR)', "\n";

print '# Load statmod, which needs to have been previously installed in the edgeR environment:', "\n";
print 'library(statmod)', "\n";

print '# Import data:', "\n";
print 'setwd("', $working_dir, '")', "\n";

# one very long line:
print 'edgeR_input_data';
print ' <- read.delim("', $input_data, '",row.names="Gene")', "\n";

# another very long line:
print 'genotypes <- factor(c(', $genotype_vector, '))', "\n";
print 'geno_comp <- DGEList(counts=edgeR_input_data,group=genotypes)', "\n";

print '# Filter out weakly expressed genes using the following commands:', "\n";
print 'keep <- filterByExpr(geno_comp)', "\n";
print 'geno_comp <- geno_comp[keep, , keep.lib.sizes=FALSE]', "\n";
print 'geno_comp$samples', "\n";

print '# Set up statistical data set:', "\n";
print 'geno_comp <- calcNormFactors(geno_comp)', "\n";
print 'geno_design <- model.matrix(~0+genotypes)', "\n";
print 'colnames(geno_design) <- levels(geno_comp$samples$group)', "\n";
print 'geno_comp <- estimateDisp(geno_comp,geno_design,robust=TRUE)', "\n";
print 'geno_comp_fit <- glmQLFit(geno_comp, geno_design, robust=TRUE)', "\n";

# now, print out many stanzas of repetitive comparisons of different genotypes:

foreach my $repetitive_command (@lines) {
    print "$repetitive_command\n";
}

# some more one-time printing
print '# Sum up:', "\n";
print 'sessionInfo()', "\n";
print '# end it:', "\n";
print 'q()', "\n";

