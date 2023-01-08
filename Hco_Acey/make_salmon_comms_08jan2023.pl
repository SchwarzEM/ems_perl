#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $infile_1 = <>) {
    chomp $infile_1;
    if ( $infile_1 =~ /\A (\S+ \/ (\S+?)) _1 \.filt1\.fq\.gz\z/xms ) {
        my $filestem = $1;
        my $desc     = $2;
        my $infile_2 = $filestem . '_2.filt1.fq.gz';
        my $output   = $desc . '_salmon';

        if (! -r $infile_1 ) {
            die "Cannot read infile 1: $infile_1\n";
        }
       	if (! -r $infile_2 ) {
       	    die	"Cannot	read infile 2: $infile_2\n";
       	}

        my $command = 'salmon --no-version-check quant '
                      . '--threads 16 --seqBias --gcBias --posBias '
                      . '--geneMap $PROJECT/Acey/2022.02.13/annots_aux/haemonchus_contortus.PRJEB506.WBPS16.cds2gene.tsv.txt '
                      . '--libType A '
                      . '--index $PROJECT/Acey/2022.02.13/Hco/salmon/dbs/Hco_WBPS16_decoys.txt_gentrome_index '
                      . "--mates1 $infile_1 "
                      . "--mates2 $infile_2 "
                      . "--output $output ;"
                      ;
        print "$command\n";
    }
    else {
        die "Cannot parse input: $infile_1\n";
    }
}

