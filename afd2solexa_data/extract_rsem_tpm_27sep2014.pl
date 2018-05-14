#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use File::Basename;

while (my $input = <>) {
    chomp $input;
    my $basename = basename $input;
    if ( $basename =~ /\A (\S+) \.genes\.results \z/xms ) { 
        my $stem = $1;
        my $outfile = $stem . '.pme_TPM.tsv.txt';
        my $err     = $stem . '.err';
        print "    cut -f 1,10 $input | perl -ne ' s/gene_id/Gene/; print; ' | wbg2fullnames.pl -n ../gene_ids/c_elegans.PRJNA13758.WS245.gene_id.tsv.txt -i - 1>$outfile 2>$err ;\n";
    }
}


