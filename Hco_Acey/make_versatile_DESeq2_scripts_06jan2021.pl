#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Getopt::Long;

use List::MoreUtils qw(uniq);

my $working_dir = q{};
my $input_data  = q{};
my $conditions  = q{};
my @genotypes   = ();
my @pairs       = ();
my $suffix      = q{};

my $help;

GetOptions ( 'working_dir=s'    => \$working_dir,
             'input_data=s'     => \$input_data,
             'conditions=s',    => \$conditions,
             'genotypes=s{,}'   => \@genotypes,
             'pairs=s{,}'       => \@pairs,
             'suffix=s'         => \$suffix,
             'help'             => \$help,   );

if ( $help or (! $working_dir) or (! $input_data ) or (! $conditions )  or (! @genotypes ) ) {
    die "Format: make_versatile_DESeq2_scripts_06jan2021.pl \n",
        "    --working_dir|-w  [working directory for R]\n",
        "    --input_data|-i   [input data file for DESeq2]\n",
        "    --conditions|-c   [conditions table for DESeq2]\n",
        "    --genotypes|-g    [list of genotypes for input data replicate columns]\n",
        "    --pairs|-p        [list of genotype pairs to be compared]\n",
        "    --suffix|-s       [suffix for output files]\n",
        "    --help|-h         [print this message]\n",
        ;
}

my @quoted_genotypes = map { q{"} . $_ . q{"} } @genotypes;
my $genotype_vector  = join ",", @quoted_genotypes;

my @lines = ();

my $i = 0;

while (@pairs) {
    my $geno1 = shift @pairs;
    my $geno2 = shift @pairs;

    my $line1 = '#' ." $geno1 vs. $geno2:";
    push @lines, $line1;

    # res_XXX_vs_YYY <- results(dds1, contrast=c("condition","XXX","YYY"), alpha = 0.1)
    # Note that I am doing it this way so that if XXX has stronger expression, it will get a positive logFC value!
    # That does not happen if one does things intuitively and submits 'pair=c("XXX","YYY")' to edgeR.

    my $line2 = 'res_' . $geno1 . '_vs_' . $geno2 . ' <- results(dds1, contrast=c("condition","' . $geno1 . '","' . $geno2 . '"), alpha = 0.1)';
    push @lines, $line2;

    # summary(res_XXX_vs_YYY)
    my $line3 = 'summary(res_' . $geno1 . '_vs_' . $geno2 . ')';
    push @lines, $line3;

    # resOrdered_XXX_vs_YYY <- res[order(res_XXX_vs_YYY$padj),]
    my $line4 = 'resOrdered_' . $geno1 . '_vs_' . $geno2 . ' <- res[order(res_' . $geno1 . '_vs_' . $geno2 . '$padj),]';
    push @lines, $line4;

    # write.table(as.data.frame(resOrdered_XXX_vs_YYY), file="XXX.vs.YYY_DESeq2_all.data_$suffix.csv", sep = "\t")
    my $line5 = 'write.table(as.data.frame(resOrdered_' 
                . $geno1 
                . '_vs_' 
                . $geno2 
                . '), file=' 
                . "\"$geno1.vs.$geno2" 
                . '_DESeq2_all.data_' 
                . $suffix 
                . ".orig.tsv.txt\""
                . ', sep="\t", col.names=NA)'
                ;
    push @lines, $line5;
}

# Print a bunch of one-time stuff:

print '# Load DESeq2:', "\n";
print 'library(DESeq2)', "\n";
print '# Import data:', "\n";
print 'setwd("', $working_dir, '")', "\n";

print 'DESeq2_input_data';
print ' <- read.delim("', $input_data, '",row.names="Gene")', "\n";

print 'DESeq2_conditions';
print ' <- read.delim("', $conditions, "\")\n";

print "dds1 <- DESeqDataSetFromMatrix(countData = DESeq2_input_data, colData = DESeq2_conditions, design = ~ batch + condition)\n";
print "dds1 <- DESeq(dds1)\n";

# now, print out many stanzas of repetitive comparisons of different genotypes:

foreach my $repetitive_command (@lines) {
    print "$repetitive_command\n";
}

# some more one-time printing

print '# Sum up:', "\n";
print 'sessionInfo()', "\n";
print '# end it:', "\n";
print 'q()', "\n";

