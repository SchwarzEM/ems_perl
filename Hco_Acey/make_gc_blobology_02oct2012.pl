#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

my $fasta_file       = 'SOAPdn.sp9.16sep2012a.fa';

my $bt2_db_file      = 'bt_dbs/SOAPdn.sp9_genDNA';
my $bt2_db_file_base = basename $bt2_db_file;

print '#!/bin/bash', "\n\n";

while (my $readfile = <>) { 
    chomp $readfile;
    my $readfile_stem = basename $readfile;
    $readfile_stem =~ s/\.jumbled\.fa//;

    print "    bowtie2 -p 12 --very-fast-local --reorder --mm -x $bt2_db_file -f -U $readfile |",
          " tee >(samtools view -S -b - > $bt2_db_file_base.$readfile_stem.bowtie2.bam) |",
          " tee >(samtools view -S -b - | samtools sort -m 2000000000 - $bt2_db_file_base.$readfile_stem.bowtie2.sorted) |",
          " sam_len_cov_gc_insert.pl -f $fasta_file -s - -out $bt2_db_file_base.$readfile_stem. ;\n";
    print "\n";
}

print "    program_done_e-ping.pl -p done_gc_blobology_01oct2012 ;\n\n";

