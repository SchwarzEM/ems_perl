#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = '#!/bin/bash' . "\n\n";

# sample input: trimmed_reads/10488_6539_90649_HLYCVBGX9_vos2num3_TTAGGC_R1.trim_filt1.fq.gz

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\Atrimmed_reads\/(\S+_R1).trim_filt1.fq.gz\z/xms ) {
        my $stem      = $1;

        my $outfile   = 'trimmed_reads_50nt/' . "$stem.50nt_trim_filt1.fq";
        my $in_count  = "$stem.trim_filt1.fq.count.txt";
        my $out_count = "$stem.50nt_trim_filt1.fq.count.txt";

        $outfile   = safename($outfile);
        $in_count  = safename($in_count);
        $out_count = safename($out_count);

        # print header only once at the start
        print $header if $header;
        $header = q{};

        print "zcat $input | quality_trim_fastq.pl -i - -u 50 -o $outfile ;\n";
        print "zcat $input | count_simple_fastq_residues.pl > $in_count ;\n";
        print "mv -i $in_count trimmed_reads ;\n";
        print "cat $outfile | count_simple_fastq_residues.pl > $out_count ;\n";
        print "mv -i $out_count trimmed_reads_50nt ;\n";
        print "gzip -1 $outfile ;\n";
        print "\n";
    }
    else {
        die "Cannot parse input line: $input\n";
    }
}


sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}


