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
        $desc =~ s/\s/./g;
        $desc =~ s/\.with\.Phusion//;

        print "    mkdir /sternlab/redivivus/data02/schwarz/Acey_genomics/$project ;\n";
        print "    cd /sternlab/redivivus/data02/schwarz/Acey_genomics/$project ;\n";

        print '    wget --user=gec --password=gecilluminadata',
              ' --output-file ../misc/', $project, '_logfile_21may2012.txt',
              ' --recursive --level=1 --no-parent --no-directories --no-check-certificate --accept .fastq.gz',
              ' https://jumpgate.caltech.edu/runfolders/volvox/120509_SN787_0116_BC0D9HACXX/Unaligned',
              '/Project_', $project, '_index', $index, '/Sample_', "$project ;\n",
              ;

        # We need to stick in a basic quality trim ("-n") to get rid of reads that didn't pass chastity!
        #    quality_trim_fastq.pl -i - -u 100 -q 33 -n -o 
        # However, *any* such filtering inexorably causes jumbling, and I want to do it *after* all the other text-wrangling, to minimize the jumble-up.
        # So...

        print '    zcat ', $project, '*_R1_*.fastq.gz > /sternlab/redivivus/data02/schwarz/Acey_genomics/raw_', $raw_type, '_reads/', $project, '_', $desc, "_orig_21may2012.R1.fq ;\n";

        print '    zcat ', $project, '*_R2_*.fastq.gz > /sternlab/redivivus/data02/schwarz/Acey_genomics/raw_', $raw_type, '_reads/', $project, '_', $desc, "_orig_21may2012.R2.fq ;\n";

        print "    cd /sternlab/redivivus/data02/schwarz/Acey_genomics ;\n";

        print '    rm -rf /sternlab/redivivus/data02/schwarz/Acey_genomics/', "$project ;\n";

        print '    shuffleSequences_fastq.pl',
              ' raw_', $raw_type, '_reads/', $project, '_', $desc, '_orig_21may2012.R1.fq',
              ' raw_', $raw_type, '_reads/', $project, '_', $desc, '_orig_21may2012.R2.fq',
              ' raw_', $raw_type, '_reads/', $project, '_', $desc, "_orig_21may2012.fq ;\n",
              ; 

        print '    rm raw_', $raw_type, '_reads/', $project, '_', $desc, '_orig_21may2012.R1.fq raw_', $raw_type, '_reads/', $project, '_', $desc, '_orig_21may2012.R2.fq', " ;\n";

        print '    retroname_fastq_reads.pl raw_', $raw_type, '_reads/', $project, '_', $desc, '_orig_21may2012.fq > raw_', $raw_type, '_reads/', $project, '_', $desc, "_retro_21may2012.fq ;\n";

        print '    rm raw_', $raw_type, '_reads/', $project, '_', $desc, "_orig_21may2012.fq ;\n";

        print '    quality_trim_fastq.pl -u 100 -q 33 -n',
              ' -i raw_', $raw_type, '_reads/', $project, '_', $desc, '_retro_21may2012.fq',
              ' -o raw_', $raw_type, '_reads/', $project, '_', $desc, "_retro_21may2012.jumbled.fq ;\n",
              ;

        print '    rm raw_', $raw_type, '_reads/', $project, '_', $desc, "_retro_21may2012.fq ;\n";

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

