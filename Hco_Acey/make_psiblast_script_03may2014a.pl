#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use File::Basename;
use Scalar::Util qw(looks_like_number);

my $e_value = $ARGV[0];

if ( (! $e_value) or ( $e_value <= 0 ) or (! looks_like_number($e_value) ) ) { 
    die "E-value needs to be positive real number, not this: $e_value\n"
}

my @queries = qw(
    queries_03may2014/Acey_s0010.g910.t3_22-171.fa 
    queries_03may2014/Acey_s0010.g910.t3_172-316.fa
    queries_03may2014/Acey_s0010.g910.t3_317-483.fa
    queries_03may2014/Acey_s0010.g910.t3_484-646.fa
    queries_03may2014/Acey_s0010.g910.t3_658-804.fa

    queries_03may2014/Acey_s0004.g1962.t1.fa
    queries_03may2014/Acey_s0004.g2039.t1.fa
    queries_03may2014/Acey_s0065.g3635.t1_184-286.fa
    queries_03may2014/Acey_s0065.g3635.t1_78-181.fa

    queries_03may2014/Acey_s0010.g996.t1.fa
    queries_03may2014/Acey_s0212.g2239.t1.fa
    queries_03may2014/Acey_s0517.g2812.t1.fa
);

my $header = '#!/bin/bash' . "\n\n";
my $footer = "\n    echo done_psiblast_script_03may2014a > script_03may2014a.txt ;\n\n";

foreach my $query (@queries) { 
    print $header if $header;
    $header = q{};
    my $stem = $query;
    $stem = basename($stem);
    $stem    =~ s/\.fa\z//;
    print "    psiblast -db dbs/metazoan_prots_14apr2014",
          " -query $query -num_threads 8 -evalue $e_value -inclusion_ethresh $e_value -num_iterations 20 -outfmt 6",
          " -seg yes -out $stem.psiblast_20x_", "$e_value.metazoan_prots.tsv.txt ;\n",
          ;
}
print $footer if (! $header);


