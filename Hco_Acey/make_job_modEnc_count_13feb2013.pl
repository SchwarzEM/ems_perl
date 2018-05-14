#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

my $header = '#!/bin/bash' 
             . "\n\n"
             . "    echo >  modEnc_count_13feb2013.txt ;\n\n"
             ;

my @fastqs = ();

while (my $input = <>) {
    chomp $input;
    push @fastqs, $input;
}

foreach my $fastq (@fastqs) { 
    my $basename = basename $fastq;

    print $header if $header;
    $header = q{};

    print "    echo >> modEnc_count_13feb2013.txt ;\n";    
    print "    echo \"    FastQ file -- $basename:\" >> modEnc_count_13feb2013.txt ;\n";
    print "    echo >> modEnc_count_13feb2013.txt ;\n";
    print "    zcat $fastq | fastq2fa_simple.pl | count_fasta_residues.pl -i - -e >> modEnc_count_13feb2013.txt ;\n";
    print "    echo >> modEnc_count_13feb2013.txt ;\n";
    print "\n";

}

if (@fastqs) {
    print "    e_ping -p done_modEnc_count_13feb2013 ;\n\n";
}


