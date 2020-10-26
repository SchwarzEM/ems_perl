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
                      . '--seqBias --gcBias --validateMappings --threads 8 '
                      . '--geneMap $SCRATCH/Acey/2020.07.28/salmon/HAEM_V4_final.tx2gene.txt '
                      . '--libType A '
                      . '--index $SCRATCH/Acey/2020.07.28/salmon/dbs/HAEM_V4_final_gentrome_index '
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


