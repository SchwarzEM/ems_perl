#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $fasta = <> ) {
    chomp $fasta;
    my $gff     = q{};
    my $cds_dna = q{};
    my $pep     = q{};
    if ( $fasta =~ /\A (\S+) \. fa \z/xms ) {
        my $stem = $1;
        $gff     = "$stem.gff";
        $cds_dna = "$stem.cds_dna.fa";
        $pep     = "$stem.pep.fa";
    }
    else {
        die "Cannot parse: $fasta\n";
    }
    $fasta = "/ocean/projects/mcb190015p/schwarze/wallacei/wal_jun2026/HKU_genomes/$fasta";
    $gff   = "/ocean/projects/mcb190015p/schwarze/wallacei/wal_jun2026/HKU_GFFs/$gff";
    if (! -e $fasta ) {
        die "Cannot find: $fasta\n";
    }
    if (! -e $gff ) {
        die "Cannot find: $gff\n";
    }
    print "gffread -g $fasta -o /dev/null --keep-genes -C -P -V -H -x $cds_dna -y $pep $gff ;\n";
}

print "mamba deactivate ;\n";

