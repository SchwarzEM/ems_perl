#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Spec::Functions;  # catdir, catfile

my $header = '#!/bin/bash';

my $old_data_dir = '/home/bioinformatics/sepal_rnaseq_dec2016/TRAP-combined/5-9-17TRAPseq/trimmed_reads';

my $new_data_dir = '/home/bioinformatics/sepal_rnaseq_dec2016/TRAP-combined/6_29_2017_TRAPseq/trimmed_reads';

my @data_types = qw( GCI1 GCI6 GCI7 GCI12 SCI4 SCI5 SCI8 );

my %old_data2file = (
    'GCI1'  => '8435_9081_52321_HJKHJBCXY_GCI1_may2017_ATCACG_R1.trim_filt1.fq.gz',
    'GCI6'  => '8435_9081_52322_HJKHJBCXY_GCI6_may2017_GCCAAT_R1.trim_filt1.fq.gz',
    'GCI7'  => '8435_9081_52323_HJKHJBCXY_GCI7_may2017_CAGATC_R1.trim_filt1.fq.gz',
    'GCI12' => '8435_9081_52324_HJKHJBCXY_GCI12_may2017_CTTGTA_R1.trim_filt1.fq.gz',
    'SCI4'  => '8435_9081_52325_HJKHJBCXY_SCI4_may2017_TGACCA_R1.trim_filt1.fq.gz',
    'SCI5'  => '8435_9081_52326_HJKHJBCXY_SCI5_may2017_ACAGTG_R1.trim_filt1.fq.gz',
    'SCI8'  => '8435_9081_52327_HJKHJBCXY_SCI8_may2017_ACTTGA_R1.trim_filt1.fq.gz',
);

my %new_data2file = (
    'GCI1'  => '8435_9081_52321_HL25LBCXY_GCI1_june2017_ATCACG_R1.trim_filt1.fq.gz',
    'GCI6'  => '8435_9081_52322_HL25LBCXY_GCI6_june2017_GCCAAT_R1.trim_filt1.fq.gz',
    'GCI7'  => '8435_9081_52323_HL25LBCXY_GCI7_june2017_CAGATC_R1.trim_filt1.fq.gz',
    'GCI12' => '8435_9081_52324_HL25LBCXY_GCI12_june2017_CTTGTA_R1.trim_filt1.fq.gz',
    'SCI4'  => '8435_9081_52325_HL25LBCXY_SCI4_june2017_TGACCA_R1.trim_filt1.fq.gz',
    'SCI5'  => '8435_9081_52326_HL25LBCXY_SCI5_june2017_ACAGTG_R1.trim_filt1.fq.gz',
    'SCI8'  => '8435_9081_52327_HL25LBCXY_SCI8_june2017_ACTTGA_R1.trim_filt1.fq.gz',
);

foreach my $data_type (@data_types) {
    my $old_infile = $old_data2file{$data_type};
    my $new_infile = $new_data2file{$data_type};

    my $test_pattern = q{_} . $data_type . q{_};
    if ( $old_infile !~ /$test_pattern/ ) { 
        die "Fail to recognize correct data type in old input file: $old_infile\n";
    }
    if ( $new_infile !~ /$test_pattern/ ) { 
        die "Fail to recognize correct data type in new input file: $new_infile\n";
    }

    my $outfile = $new_infile;
    $outfile =~ s/_H\w+Y_/_/;
    $outfile =~ s/june2017/may.and.june2017/;
    # remember that the output will *not* be gzipped at first
    $outfile =~ s/\.gz//;
    if (-e $outfile ) {
        die "There already exists a file with the planned output file name: $outfile\n";
    }

    my $readcount_file = "$outfile.count.txt";

    $old_infile = catfile($old_data_dir, $old_infile);
    $new_infile = catfile($new_data_dir, $new_infile);

    if (! -r $old_infile ) {
        die "Cannot read putative older input file: $old_infile\n";
    }

    if (! -r $new_infile ) {
        die "Cannot read putative newer input file: $new_infile\n";
    }

    print "$header\n\n" if $header;
    $header = q{};

    print "    zcat $old_infile $new_infile > $outfile ;\n";
    print "    cat $outfile | count_simple_fastq_residues.pl > $readcount_file ;\n";
    print "    gzip -1 $outfile ;\n";
    print "\n";
}

