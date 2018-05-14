#!/usr/bin/env perl

# make_retro_reads_14nov2013.pl -- Erich Schwarz <ems394@cornell.edu>, 11/14/2013.

use strict;
use warnings;

my $work_dir  = '/woldlab/loxcyc/data00/schwarz/RNAseq_data.16';
my $work_date = '14nov2013';

my $data_ref;

# Damn the torpodoes.
my $cpus = 12;

my @libraries = ( 13767, 13768, 13769, 13770 ) ;

# 13767) Index #9 CEM_DL_61113
# 13768) Index #10 CEM_DR_61113
# 13769) Index #11 CEM_VL_61113
# 13770) Index #12 CEM_VR_61113

$data_ref->{'lib'}->{13767}->{'index'} = 9;
$data_ref->{'lib'}->{13768}->{'index'} = 10;
$data_ref->{'lib'}->{13769}->{'index'} = 11;
$data_ref->{'lib'}->{13770}->{'index'} = 12;

$data_ref->{'lib'}->{13767}->{'flowcell'} = 'C2728ACXX'; 
$data_ref->{'lib'}->{13768}->{'flowcell'} = 'C2728ACXX';
$data_ref->{'lib'}->{13769}->{'flowcell'} = 'C2728ACXX';
$data_ref->{'lib'}->{13770}->{'flowcell'} = 'C2728ACXX';

$data_ref->{'lib'}->{13767}->{'description'} = 'CEM_DL_pool';
$data_ref->{'lib'}->{13768}->{'description'} = 'CEM_DR_pool';
$data_ref->{'lib'}->{13769}->{'description'} = 'CEM_VL_pool';
$data_ref->{'lib'}->{13770}->{'description'} = 'CEM_VR_pool';

$data_ref->{'lib'}->{13767}->{'data_URL'} = 'https://jumpgate.caltech.edu/runfolders/volvox/130620_SN787_0148_B';
$data_ref->{'lib'}->{13768}->{'data_URL'} = 'https://jumpgate.caltech.edu/runfolders/volvox/130620_SN787_0148_B';
$data_ref->{'lib'}->{13769}->{'data_URL'} = 'https://jumpgate.caltech.edu/runfolders/volvox/130620_SN787_0148_B';
$data_ref->{'lib'}->{13770}->{'data_URL'} = 'https://jumpgate.caltech.edu/runfolders/volvox/130620_SN787_0148_B';

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
          " $data_URL$flowcell", '/Unaligned.I6/Project_', $library, '_index', $index, '/Sample_', "$library ;\n",
          ;

    # We need to stick in a basic quality trim ("-n") to get rid of reads that didn't pass chastity!
    # For single-end data, no jumbling.  Yay!
    # But I do want to have both 38-nt reads (exactly) and 50-nt reads (exactly).

    print '    zcat ', $library, '*_R1_*.fastq.gz > ', $work_dir, '/', $desc, '_orig_', "$work_date.R1.fq ;\n";

    print "    cd $work_dir ;\n";
    print "    rm -rf $work_dir/$library ;\n";

    print '    retroname_fastq_reads.pl ', $desc, '_orig_', "$work_date.R1.fq > $desc", '_pre.retro_', "$work_date.R1.fq ;\n";
    print "    rm $desc", '_orig_', "$work_date.R1.fq ;\n";

    print "    quality_trim_fastq.pl -q 33 -u 50 -n -t 3 -m 50 -i $desc", '_pre.retro_', "$work_date.R1.fq -o $desc", '_retro50nt_', "$work_date.R1.fq ;\n";
    print "    quality_trim_fastq.pl -q 33 -u 38 -n -t 3 -m 38 -i $desc", '_pre.retro_', "$work_date.R1.fq -o $desc", '_retro38nt_', "$work_date.R1.fq ;\n";

    print '    e_ping -p got_initial_', "$desc.RNAseq_reads.$work_date ;\n";

    print "\n";
}

foreach my $library (@libraries) {
    my $index    = $data_ref->{'lib'}->{$library}->{'index'};
    my $flowcell = $data_ref->{'lib'}->{$library}->{'flowcell'};
    my $desc     = $data_ref->{'lib'}->{$library}->{'description'};
    my $data_URL = $data_ref->{'lib'}->{$library}->{'data_URL'};

    print "    bowtie ../ce6splices/ce6spl38spk -p $cpus -k 11 -v 2 --best -m 10 -t --strata",
          " -q $desc", '_retro38nt_', "$work_date.R1.fq $desc", '_retro38nt_', "$work_date.ce6spl38spk.bowtie.txt ;\n",
          ;

    print '    e_ping -p done_', "$desc", '_retro38nt_', "$work_date.ce6spl38spk.bowtie ;\n";
    print "\n";

    print '     python $ERANGEPATH/makerdsfrombowtie.py', " $desc", '_retro38nt_', "$work_date",
          " $desc", '_retro38nt_', "$work_date.ce6spl38spk.bowtie.txt",
          " $desc", '_retro38nt_', "$work_date.ce6spl38spk.bt.rds -RNA ../ce6splices/sangerGene.txt -cache 1000000 -index ;\n",
          ;

    print '    e_ping -p done_', "$desc", '_retro38nt_', "$work_date.ce6spl38spk.bt.rds ;\n";
    print "\n";

    print "    mv $desc", '_retro38nt_', "$work_date.ce6spl38spk.bt.wig $desc", '_retro38nt_', "$work_date.ce6spl38spk.bt.wig.orig ;\n";
    print "    censor_worm_wigs.pl $desc", '_retro38nt_', "$work_date.ce6spl38spk.bt.wig.orig > $desc", '_retro38nt_', "$work_date.ce6spl38spk.bt.wig ;\n";
    print '    e_ping -p done_', "$desc", '_retro38nt_', "$work_date.makewiggle ;\n";
    print "\n";

    print "    ../ce6splices/runWormAnalysis_ems.sh $desc", '_retro38nt_', "$work_date.ce6spl38spk.bt ;\n";
    print "    mv -i $desc", '_retro38nt_', "$work_date.ce6spl38spk.bt.final.rpkm $desc", '_retro38nt_', "$work_date.ce6spl38spk.bt.final.WS190.rpkm ;\n";
    print '    update_wbg_ids.pl -t /sternlab/redivivus/data02/schwarz/archive/taygeta/science/afd2solexa_paper/ws220_xace_data/ws220_old2new_wbgenes.tsv',
          " $desc", '_retro38nt_', "$work_date.ce6spl38spk.bt.final.WS190.rpkm > $desc", '_retro38nt_', "$work_date.ce6spl38spk.bt.final.WS220.rpkm ;\n";
    print '    part2wbfullgenes.pl --names /sternlab/redivivus/data02/schwarz/archive/taygeta/science/afd2solexa_paper/analysis_11nov2010/WS220_WBGeneIDs.tsv',
          " $desc", '_retro38nt_', "$work_date.ce6spl38spk.bt.final.WS220.rpkm 1>",
          "$desc", '_retro38nt_', "$work_date.ce6spl38spk.bt.WS220.rpkm 2>",
          "$desc", '_retro38nt_', "$work_date.ce6spl38spk.bt.WS220.err.txt ;\n",
          ;
    print "    rm $desc", '_retro38nt_', "$work_date.ce6spl38spk.bt.final.WS220.rpkm ;\n";
    print '    e_ping -p done_', "$desc", '_retro38nt_', "$work_date.runWormAnalysis ;\n";
    print "\n";
}

