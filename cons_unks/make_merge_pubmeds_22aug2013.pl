#!/usr/bin/env perl

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A prot2pubmeds\/ (\S+) _prot2pubmed\.alleged\.txt \z/xms ) { 
        my $prefix = $1;
        print '    cat prot2pubmeds/', $prefix, '_prot2pubmed.alleged.txt prot2pubmeds/', $prefix, '_prot2pubmed.direct.txt',
              ' | /mnt/home/emsch/perl.svn/trunk/cons_unks/merge_alt_prot2pubmeds.pl 1>prot2pubmeds/', $prefix, '_prot2pubmed.merged.txt',
              ' 2>', $prefix, "_prot2pubmed.merged.err ;\n",
              ;
    }
}

