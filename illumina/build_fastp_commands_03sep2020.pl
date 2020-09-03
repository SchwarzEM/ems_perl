#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input_1 = <>) {
    chomp $input_1 ;
    if ( $input_1 =~ /\A (\S+  \/ ([^\s\/]+)) _1\.fq\.gz \z/xms ) {
        my $full_stem = $1;
        my $file_stem = $2;

        my $input_2  = $full_stem . '_2.fq.gz';

        my $output_1 = $file_stem . '_1.filt1.fq.gz';
       	my $output_2 = $file_stem . '_2.filt1.fq.gz';
        my $unpaired = $file_stem . '.unpaired.fq.gz';

        my $output_1_count = $file_stem . '_1.filt1.fq.count.txt';
        my $output_2_count = $file_stem . '_2.filt1.fq.count.txt';
        my $unpaired_count = $file_stem . '.unpaired.fq.count.txt';

        my $json = $file_stem . '.json';
        my $html = $file_stem . '.html';

        if (! -r $input_1 ) {
            die "Can't read input 1: $input_1\n";
        }
        if (! -r $input_2 ) {
            die "Can't read input 2: $input_2\n";
        }

        my $command1 = "fastp --thread 8 --dont_overwrite --length_required 125 --max_len1 150 "
                       . "--json $json --html $html "
                       . "--in1 $input_1 --out1 $output_1 "
                       . "--in2 $input_2 --out2 $output_2 "
                       . "-unpaired1 $unpaired --unpaired2 $unpaired ;"
                       ;

        my $command2 = "zcat $output_1 | count_simple_fastq_residues.pl > $output_1_count ;";
       	my $command3 = "zcat $output_2 | count_simple_fastq_residues.pl > $output_2_count ;";
       	my $command4 = "zcat $unpaired | count_simple_fastq_residues.pl > $unpaired_count ;";

        print "$command1\n";
        print "$command2\n"; 
        print "$command3\n"; 
        print "$command4\n"; 
        print "\n";
    }
    else {
        die "Cannot parse: $input_1\n";
    }
}

