#!/usr/bin/env perl

use strict;
use warnings;

my @conditions = (
    [ ("ACEY.L3i", "ACEY.24HCM") ],
    [ ("ACEY.L3i", "ACEY.24.PI") ],
    [ ("ACEY.24HCM", "ACEY.24.PI") ],
    [ ("ACEY.24.PI", "ACEY.5.D") ],
    [ ("ACEY.5.D", "ACEY.post_12.D") ],
    [ ("ACEY_18D.noAlb.4hr", "ACEY_18D.Alb.4hr") ],
    [ ("ACEY_18D.HEPES.4hr", "ACEY_18D.Cry5B.4hr") ],
    [ ("ACEY_18D.HEPES.24hr", "ACEY_18D.Cry5B.24hr") ],
    [ ("ACEY_18D.SB.plusHEPES", "ACEY_18D.SB.plusCry5B") ],
    [ ("ACEY_18D.HEPES.24hr", "ACEY_18D.SB.plusHEPES") ],
    [ ("ACEY_18D.Cry5B.24hr", "ACEY_18D.SB.plusCry5B") ],
);

my $cond_count = @conditions;
$cond_count--;

print "library(DESeq)\n";
print "setwd(\"/sternlab/redivivus/data02/schwarz/Acey_genomics/ng_revision_mar2014/rnaseq\")\n";
print "AceyCountTable = read.delim(\"Acey_pme_expected_count.rounded.5plus_reads.28mar2014.txt\",row.names=\"Gene\")\n";

# For this version, we are putting ACEY.12.D, ACEY.17.D and ACEY.19.D into a single category: ACEY.post_12.D.
# As before, We use '_' instead of '-' in labels to keep R from choking on them, but the real problem was probably failing to put in commas!
print "condition = factor( c(",
      " \"ACEY.L3i\", \"ACEY.24HCM\", \"ACEY.24.PI\", \"ACEY.5.D\", \"ACEY.post_12.D\", \"ACEY.post_12.D\",", 
      " \"ACEY_18D.Alb.4hr\", \"ACEY_18D.noAlb.4hr\", \"ACEY_18D.Cry5B.4hr\", \"ACEY_18D.HEPES.4hr\", \"ACEY_18D.SB.plusCry5B\",",
      " \"ACEY_18D.SB.plusHEPES\", \"ACEY_18D.Cry5B.24hr\", \"ACEY_18D.HEPES.24hr\", \"ACEY.post_12.D\" ) )\n",
      ;

print "cds = newCountDataSet( AceyCountTable, condition )\n";
print "cds = estimateSizeFactors( cds )\n";
print "cds = estimateDispersions( cds, method=\"pooled\", sharingMode=\"maximum\" )\n";

foreach my $index (0..$cond_count--) { 
    my $i = $index;
    $i++;
    my @conds = @{ $conditions[$index] };
    my $cond_a = $conds[0];
    my $cond_b = $conds[1];
    print "res$i = nbinomTest( cds, \"$cond_a\", \"$cond_b\" )\n";
    print "write.table(res$i, file=\"DESeq.comp.$cond_a.vs.$cond_b.30mar2013.post_12.D_trip.raw.txt\")\n";
}
print "q()\n";

