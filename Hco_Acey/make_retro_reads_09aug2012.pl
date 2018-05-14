#!/usr/bin/env perl

# make_retro_reads_21may2012.pl -- Erich Schwarz <emsch@caltech.edu>, 5/21/2012.

use strict;
use warnings;
use Getopt::Long;

my $infile = q{};
my $raw_type = q{};
my $bash;

GetOptions ( 'infile=s'    => \$infile,
             'raw_type=s', => \$raw_type, 
             'bash'        => \$bash, );


if (! -e '/sternlab/redivivus/data02/schwarz/Acey_genomics/misc') {
    die "Can't find required directory /sternlab/redivivus/data02/schwarz/Acey_genomics/misc\n";
}

if ( $raw_type !~ /\A \w+ \z/xms ) { 
    die "Can't parse raw_type \"$raw_type\"\n";
}

open my $INFILE, '<', $infile or die "Can't open infile $infile: $!";

if ($bash) { 
    print '#!/bin/bash', "\n";
}

print "\n";

while (my $input = <$INFILE>) { 
    chomp $input;

    # Sample inputs:
    # Lane #1 : (12665) Index #9 ACEY 12 D with Phusion (None)
    # Lane #2 : (12723) Index #3 [A. ceylanicum] L3i (None)
    # Lane #2 : (12805) Index #6 ACEY 4 JUMP 2 (None)

    if ( $input =~ /\A Lane [ ] [#] \d+ [ ] : [ ] \( (\d+) \) [ ] Index [ ] [#] (\d+) [ ] (\S.+\S) [ ] \( None \) \s* \z/xms ) { 
        my $project = $1;
        my $index   = $2;
        my $desc    = $3;

        $desc =~ s/\[A\.[ ]ceylanicum\]/ACEY/g;
        $desc =~ s/Acey/ACEY/g;
        $desc =~ s/\s/./g;
        $desc =~ s/\.with\.Phusion//;

        print "    mkdir /sternlab/redivivus/data02/schwarz/Acey_genomics/$project ;\n";
        print "    cd /sternlab/redivivus/data02/schwarz/Acey_genomics/$project ;\n";

        print '    wget --user=gec --password=gecilluminadata',
              ' --output-file ../misc/', $project, '_logfile_21may2012.txt',
              ' --recursive --level=1 --no-parent --no-directories --no-check-certificate --accept .fastq.gz',
              ' https://jumpgate.caltech.edu/runfolders/Volvox/120803_SN787_0127_BD17FPACXX/Unaligned',
              '/Project_', $project, '_index', $index, '/Sample_', "$project ;\n",
              ;

        # We need to stick in a basic quality trim ("-n") to get rid of reads that didn't pass chastity!
        # I also want to enforce a 40-nt minimum on reads.
        # Note minor point: for genomic assembly I enforced '-t 3' for quality, but for the RNA-seq reads I didn't bother.
        # Moreover, I only remembered that issue *now*, so, for the sake of consistency, don't bother this time either.

        #    quality_trim_fastq.pl -q 33 -u 50 -m 40 -n

        # For single-end data, no jumbling.  Yay!
        # So...

        print '    zcat ', $project, '*_R1_*.fastq.gz > /sternlab/redivivus/data02/schwarz/Acey_genomics/raw_', $raw_type, '_reads/', $project, '_', $desc, "_orig_06aug2012.R1.fq ;\n";

        # No need to do this for single-read data:
        # print '    zcat ', $project, '*_R2_*.fastq.gz > /sternlab/redivivus/data02/schwarz/Acey_genomics/raw_', $raw_type, '_reads/', $project, '_', $desc, "_orig_06aug2012.R2.fq ;\n";

        print "    cd /sternlab/redivivus/data02/schwarz/Acey_genomics ;\n";

        print '    rm -rf /sternlab/redivivus/data02/schwarz/Acey_genomics/', "$project ;\n";

        # No need to do this for single-read data:
        # print '    shuffleSequences_fastq.pl',
        #       ' raw_', $raw_type, '_reads/', $project, '_', $desc, '_orig_21may2012.R1.fq',
        #       ' raw_', $raw_type, '_reads/', $project, '_', $desc, '_orig_21may2012.R2.fq',
        #      ' raw_', $raw_type, '_reads/', $project, '_', $desc, "_orig_21may2012.fq ;\n",
        #       ; 
        # 
        # print '    rm raw_', $raw_type, '_reads/', $project, '_', $desc, '_orig_21may2012.R1.fq raw_', $raw_type, '_reads/', $project, '_', $desc, '_orig_21may2012.R2.fq', " ;\n";

        print '    retroname_fastq_reads.pl raw_', $raw_type, '_reads/', $project, '_', $desc, '_orig_06aug2012.R1.fq',
              ' > ', 
              'raw_', $raw_type, '_reads/', $project, '_', $desc, "_pre.retro_06aug2012.R1.fq ;\n";

        print '    rm raw_', $raw_type, '_reads/', $project, '_', $desc, "_orig_06aug2012.R1.fq ;\n";

        print '    quality_trim_fastq.pl -q 33 -u 50 -m 40 -n',
              ' -i raw_', $raw_type, '_reads/', $project, '_', $desc, '_pre.retro_06aug2012.R1.fq',
              ' -o raw_', $raw_type, '_reads/', $project, '_', $desc, "_retro_06aug2012.R1.fq ;\n",
              ;

        print '    rm raw_', $raw_type, '_reads/', $project, '_', $desc, "_pre.retro_06aug2012.R1.fq ;\n";

        print '    bowtie filter_seqs/CriGri1_AceyMito_IlluPrim -p 6 -t',
              ' --un filt2_', $raw_type, '_reads_23jul2012/', $project, '_', $desc, '_retro_06aug2012.no_cont.min40.R1.se.fq',
              ' -q raw_', $raw_type, '_reads/', $project, '_', $desc, '_retro_06aug2012.R1.fq',
              " /dev/null ;\n",
              ;

        print '    rm raw_', $raw_type, '_reads/', $project, '_', $desc, "_retro_06aug2012.R1.fq ;\n";

        print '    program_done_e-ping.pl -p got_initial_', "$project ;\n";

        print "\n";
    }
    else { 
        if ( $input =~ /\A Lane /xms ) {
            die "Can't parse input line: $input\n";
        }
    }
}

close $INFILE or die "Can't close filehandle to infile $infile: $!";

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

