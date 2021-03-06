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
    print "tbl2asn";
    print " -p .";
    print " -t /mnt/ls15/scratch/users/emsch/2017/caenogens/genbank/genome/nigoni_genbank_annot_files/nigoni_genome_template.sbt";
    print " -i $ncbi_fasta";

    # There is no point in invoking a non-existent TBL file!  Small scaffolds/contigs may well lack a TBL.
    if (-e $ncbi_table) {
        print " -f $ncbi_table";
    }

    print ' -j "[organism=Caenorhabditis nigoni] [tech=wgs] [strain=JU1422]"';
    print " -a r1u -c f -V vb -M n -X E -l align-genus -l align-trnscpt";
    print " -Z discrep_", $basename, ".txt ;";
    print "\n";
}

print "\n";

