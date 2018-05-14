#!/usr/bin/env perl

use strict;
use warnings;

my @conditions = (
    [ ("ACEY.L3i", "ACEY.24HCM") ],
    [ ("ACEY.L3i", "ACEY.24.PI") ],
    [ ("ACEY.24HCM", "ACEY.24.PI") ],
    [ ("ACEY.24.PI", "ACEY.5.D") ],
    [ ("ACEY.5.D", "ACEY.12.D") ],
    [ ("ACEY.12.D", "ACEY.17.D") ],
    [ ("ACEY.17.D", "ACEY.19.D") ],
    [ ("ACEY.18D.noAlb.4hr", "ACEY.18D.Alb.4hr") ],
    [ ("ACEY.18D.HEPES.4hr", "ACEY.18D.Cry5B.4hr") ],
    [ ("ACEY.18D.HEPES.24hr", "ACEY.18D.Cry5B.24hr") ],
    [ ("ACEY.18D.SB.plusHEPES", "ACEY.18D.SB.plusCry5B") ],
    [ ("ACEY.18D.HEPES.24hr", "ACEY.18D.SB.plusHEPES") ],
    [ ("ACEY.18D.Cry5B.24hr", "ACEY.18D.SB.plusCry5B") ],
);

my $cond_count = @conditions;
$cond_count--;

print "library(NOISeq)\n";
print "setwd(\"/mnt/home/emsch/work/Acey/ng_revision_mar2014/rnaseq\")\n";
print "AceyCountTable = read.delim(\"Acey_pme_expected_count.rounded.5plus_reads.28mar2014.txt\",row.names=\"Gene\")\n";

print "AceyFactors = data.frame(Acey_condition = c( \"ACEY.L3i\", \"ACEY.24HCM\", \"ACEY.24.PI\", \"ACEY.5.D\", \"ACEY.12.D\", \"ACEY.17.D\",",
      " \"ACEY.18D.Alb.4hr\", \"ACEY.18D.noAlb.4hr\", \"ACEY.18D.Cry5B.4hr\", \"ACEY.18D.HEPES.4hr\", \"ACEY.18D.SB.plusCry5B\", \"ACEY.18D.SB.plusHEPES\",",
      " \"ACEY.18D.Cry5B.24hr\", \"ACEY.18D.HEPES.24hr\", \"ACEY.19.D\" ) )\n",
      ;

print "AceyData <- readData(data = AceyCountTable, factors = AceyFactors)\n";

foreach my $index (0..$cond_count--) {
    my @conds = @{ $conditions[$index] };
    my $cond_a = $conds[0];
    my $cond_b = $conds[1];

    # Note that we have to list condition B *first* to get NOISeq to call 'up' genes those that go up into condition B.

    print "AceyNoiseq.$cond_a.to.$cond_b",
          " <- noiseq(AceyData, k = 0.5, norm=\"tmm\", replicates=\"no\", factor=\"Acey_condition\",",
          " conditions=c( \"$cond_b\", \"$cond_a\"), pnr = 0.2, nss = 5, v = 0.02, lc = 0)\n",
          ;

    print "AceyNoiseq.$cond_a.to.$cond_b.deg.up   = degenes(AceyNoiseq.$cond_a.to.$cond_b, q = 0.99, M = \"up\")\n";
    print "AceyNoiseq.$cond_a.to.$cond_b.deg.down = degenes(AceyNoiseq.$cond_a.to.$cond_b, q = 0.99, M = \"down\")\n";

    print "write.table(AceyNoiseq.$cond_a.to.$cond_b.deg.up, file=\"AceyNoiseq.$cond_a.to.$cond_b.deg.up.raw.txt\")\n";
    print "write.table(AceyNoiseq.$cond_a.to.$cond_b.deg.down, file=\"AceyNoiseq.$cond_a.to.$cond_b.deg.down.raw.txt\")\n";
}
print "q()\n";

