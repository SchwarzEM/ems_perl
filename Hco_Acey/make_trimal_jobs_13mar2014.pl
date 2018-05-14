#!/usr/bin/env perl

use strict;
use warnings;

my $header = '#!/bin/bash' . "\n\n";

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A muscle_aligns \/ ( \S+ _13mar2014 \.musc_) orig\.fa \z/xms ) {
        my $stem = $1;
        my $trimal_out1  = 'muscle_aligns/'  . $stem . 'gt0.75_rs0.75.fa';
        my $trimal_html1 = 'trimal_outputs/' . $stem . 'gt0.75_rs0.75.html';
        my $trimal_out2  = 'muscle_aligns/'  . $stem . 'gappyout.rs0.75.fa';
        my $trimal_html2 = 'trimal_outputs/' . $stem . 'gappyout.rs0.75.html';
        print $header if $header;
        $header = q{};
        print "    trimal -in $input -out $trimal_out1 -gt 0.75 -resoverlap 0.75 -seqoverlap 75 -keepheader -htmlout $trimal_html1 ;\n";
        print "    trimal -in $input -out $trimal_out2 -gappyout -resoverlap 0.75 -seqoverlap 75 -keepheader -htmlout $trimal_html2 ;\n";
        print "\n";
    }
    else { 
        die "Can't parse input: $input\n";
    }
}

print "    e_ping -p trimal_jobs_13mar2014 ;\n\n";


