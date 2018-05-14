#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

my $header = '#!/bin/bash' . "\n\n";

while (my $ncbi_fasta = <>) { 
    chomp $ncbi_fasta;
    my $basename = basename($ncbi_fasta);

    my $ncbi_table = $ncbi_fasta;
    $ncbi_table    =~ s/\.fsa\z/.tbl/;

    print $header if $header;
    $header = q{};
    print "tbl2asn -p . -t /mnt/ls15/scratch/users/emsch/2017/caenogens/genbank/genome/nigoni_genbank_annot_files/nigoni_genome_template.sbt",
          " -i $ncbi_fasta -f $ncbi_table",
          ' -j "[organism=Caenorhabditis nigoni] [tech=wgs] [strain=JU1422]"',
          " -a r1u -c f -V vb -M n -X E",
          " -Z discrep_", $basename, ".txt ;",
          "\n",
          ;
}

print "\n";

