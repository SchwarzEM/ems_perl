#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

my $header = "library(stats)";

# Sample input lines:
# Comparison	DESeq2_genes	edgeR_genes	Overlap_genes	Overlap_ratio
# atml1-3.vs.lgo-2.down	13	0	0	n/a
# atml1-3.vs.lgo-2.up	1	0	0	n/a
# ATML1__LGO_atml1-3.vs.atml1-3.down	89	25	23	0.920 [23/25]

while (my $input = <>) {
    chomp $input;
    # Process all lines *except* the header line.
    if ( ( $input =~ /\A (\S+) \t (\S+) \t (\S+) \t (\S+) \t [^\t]* \z/xms ) and ( $1 ne 'Comparison' ) ) { 
        my $comparison         = $1;
        my $deseq2_gene_count  = $2;
        my $edgeR_gene_count   = $3;
        my $overlap_gene_count = $4;

        # Strip out any commas place into the number for readability, then verify that they're numbers:
        $deseq2_gene_count  =~ s/[,]//g;
        $edgeR_gene_count   =~ s/[,]//g;
        $overlap_gene_count =~ s/[,]//g;

        if (    (! looks_like_number($deseq2_gene_count)  ) 
             or (! looks_like_number($edgeR_gene_count)   )
             or (! looks_like_number($overlap_gene_count) ) ) {
            die "Cannot parse one or more gene counts in: $input\n";
        }

        # get the genomewide frequency of *edgeR*-positive genes:
        my $genomewide_edgeR = ($edgeR_gene_count/27416);

        # Print the header only once, at the start:
        print "$header\n" if $header;
        $header = q{};

        # Go ahead and print this stanza even if $overlap_gene_count == 0, because this is a *two-tailed* test!
        print "# Comparison: \"$comparison\"; p-value for freq. of edgeR overlaps in the DESeq2 gene set\n";
        print 'binom.test('
            . $overlap_gene_count
            . q{,}
            . $deseq2_gene_count
            . ',p='
            . $genomewide_edgeR
            . ',alternative=c("two.sided"),conf.level=0.99)'
            . "\n"
            ;

    }
    elsif ( $input !~ /\A Comparison \t DESeq2_genes \t edgeR_genes \t Overlap_genes \t Overlap_ratio \z/xms ) { 
        die "Cannot parse input line: $input\n";
    }
}

# Print the final line only once, *if* we have actually printed data:
print "q()\n" if (! $header);

