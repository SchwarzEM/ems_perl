#!/usr/bin/env perl

# make_get_new_pristine_and_filt1_CEM_reads_27feb2017.pl -- Erich Schwarz <ems394@cornell.edu>, 2/27/2017.

use strict;
use warnings;
use autodie;

my $work_dir  = '/sternberg/gluster/gluster03/pub/schwarz/RNAseq_data.17';
my $work_date = '21feb2017';

my $data_ref;

my $cpus = 12;

my @libraries = ( 13770 ) ;

# 13770) Index #12 CEM_VR_61113

$data_ref->{'lib'}->{13770}->{'index'} = 12;

$data_ref->{'lib'}->{13770}->{'flowcell'} = 'H7T32BCXY';

$data_ref->{'lib'}->{13770}->{'description'} = 'CEM_VR_pool';

$data_ref->{'lib'}->{13770}->{'data_URL'} = 'https://jumpgate.caltech.edu/runfolders/volvox02/170224_SN787_0636_A';

print '#!/bin/bash', "\n";

print "\n";

foreach my $library (@libraries) { 
    my $index    = $data_ref->{'lib'}->{$library}->{'index'};
    my $flowcell = $data_ref->{'lib'}->{$library}->{'flowcell'};
    my $desc     = $data_ref->{'lib'}->{$library}->{'description'};
    my $data_URL = $data_ref->{'lib'}->{$library}->{'data_URL'};

    print "    mkdir $work_dir/$library ;\n";
    print "    cd $work_dir/$library ;\n";

    print '    wget --user=gec --password=gecilluminadata --output-file ../', $library, '_logfile_', $work_date, '.txt',
          ' --recursive --level=1 --no-parent --no-directories --no-check-certificate --accept .fastq.gz',
          " $data_URL$flowcell", '/Unaligned.singleIndex/Project_', $library, '_index', $index, '/Sample_', "$library ;\n",
          ;

    # We need to stick in a basic quality trim ("-n") to get rid of reads that didn't pass chastity!
    # For single-end data, no jumbling.  Yay!
    # The input reads are 100-nt, raw.  So *first* make a set trimmed for everything but max. length 50 nt, *then* make a set with
    #    max. len. 50-nt imposed as a final requirement.


    my $raw_output_file      = $desc . '_orig_' . "$work_date.R1.fq";
    my $prefilt1_output_file = $desc . '_orig_' . "$work_date.R1.len_50to100_filt1.fq";
    my $filt1_output_file    = $desc . '_orig_' . "$work_date.R1.filt1.fq";

    print '    zcat ', $library, '*_R1_*.fastq.gz > ', $work_dir, '/', "$raw_output_file ;\n";

    print "    cd $work_dir ;\n";
    print "    rm -rf $work_dir/$library ;\n";

    print "    quality_trim_fastq.pl -q 33 -n -t 3 -m 50 -i $raw_output_file      -o $prefilt1_output_file ;\n";
    print "    quality_trim_fastq.pl -q 33 -u 50         -i $prefilt1_output_file -o $filt1_output_file ;\n";

    print "\n\n";
}

