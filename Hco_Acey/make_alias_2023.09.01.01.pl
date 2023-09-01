#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;
use File::Spec::Functions;  # catdir

my $start_dir = getcwd;

while (my $input = <> ) {
    chomp $input;

    if (! -e $input ) {
        die "Apparently, this input file does not exist: $input\n";
    }

    # Sample input:
    # /ocean/projects/mcb190015p/fergusoa/Herbert_Lab_Nb_RNAseq_data/AF37_STAT6G6_WTG7_fastq/01_WT_G7_trm.fastq.gz
    if ( $input =~ /\A \/ \S+ \/ \d+ _ (\S+) _trm\.fastq\.gz /xms ) { 
        my $sub_dir  = $1;

        # This is completely arbitrary, but it's necessary to deal with an inconsistently named replicate set
        #    that is now locked into an official BioProject name.
        if ( $sub_dir eq 'STAT6_G6' ) {
            $sub_dir = 'STAT6KO_G6';
        }

        my $work_dir = catdir($start_dir, $sub_dir);

        if (! -e $work_dir ) {
            die "Apparently, this intended work directory does not exist: $work_dir\n";
        }

        print "cd $work_dir ;\n";
        print "ln -s $input ;\n";
    }
    else {
        die "Cannot parse input: $input\n";
    }
}

