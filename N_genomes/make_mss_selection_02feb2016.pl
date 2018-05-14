#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \s .* \s (\S+)\z/xms ) { 
        my $gene_name = $1;
        my $seq_name  = $2;
        print "    extract_fasta_subset.pl -f sets_of_seqs/caeno_mss_mrp_candidates_13jan2016.fa -l $seq_name";
        print " | perl -ne ' s\/$seq_name\/$gene_name  $seq_name\/; print; ' >> Cni_Csin_Crem_Cbren_MSS_seqs_02feb2016.fa ;\n";
    }
}


