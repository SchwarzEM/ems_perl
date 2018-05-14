#!/usr/bin/env perl

# make_get_new_pristine_and_filt1_CEM_reads_21feb2017.pl -- Erich Schwarz <ems394@cornell.edu>, 2/21/2017.

use strict;
use warnings;
use autodie;

my $work_dir  = '/sternberg/gluster/gluster03/pub/schwarz/RNAseq_data.17';
my $work_date = '21feb2017';

my $data_ref;

# Damn the torpedoes.
my $cpus = 12;

my @libraries = ( 13767, 13768, 13769 ) ;

# 13767) Index #9 CEM_DL_61113
# 13768) Index #10 CEM_DR_61113
# 13769) Index #11 CEM_VL_61113


$data_ref->{'lib'}->{13767}->{'index'} = 9;
$data_ref->{'lib'}->{13768}->{'index'} = 10;
$data_ref->{'lib'}->{13769}->{'index'} = 11;

$data_ref->{'lib'}->{13767}->{'flowcell'} = 'H772TBCXY'; 
$data_ref->{'lib'}->{13768}->{'flowcell'} = 'H772TBCXY';
$data_ref->{'lib'}->{13769}->{'flowcell'} = 'H772TBCXY';

$data_ref->{'lib'}->{13767}->{'description'} = 'CEM_DL_pool';
$data_ref->{'lib'}->{13768}->{'description'} = 'CEM_DR_pool';
$data_ref->{'lib'}->{13769}->{'description'} = 'CEM_VL_pool';

$data_ref->{'lib'}->{13767}->{'data_URL'} = 'https://jumpgate.caltech.edu/runfolders/volvox02/170216_SN787_0632_B';
$data_ref->{'lib'}->{13768}->{'data_URL'} = 'https://jumpgate.caltech.edu/runfolders/volvox02/170216_SN787_0632_B';
$data_ref->{'lib'}->{13769}->{'data_URL'} = 'https://jumpgate.caltech.edu/runfolders/volvox02/170216_SN787_0632_B';

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
    # But I do want to have both 38-nt reads (exactly) and 50-nt reads (exactly).

    my $raw_output_file = $desc . '_orig_' . "$work_date.R1.fq";
    my $filt1_output_file = $desc . '_orig_' . "$work_date.R1.filt1.fq";

    print '    zcat ', $library, '*_R1_*.fastq.gz > ', $work_dir, '/', "$raw_output_file ;\n";

    print "    cd $work_dir ;\n";
    print "    rm -rf $work_dir/$library ;\n";

    print "    quality_trim_fastq.pl -q 33 -n -u 50 -t 3 -m 50 -i $raw_output_file -o $filt1_output_file ;\n";

    print "\n\n";
}

