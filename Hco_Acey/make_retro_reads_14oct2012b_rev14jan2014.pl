#!/usr/bin/env perl

# make_retro_reads_14oct2012b_rev14jan2014.pl -- Erich Schwarz <ems394@cornell.edu>, 1/14/2014.

use strict;
use warnings;
use Getopt::Long;

my $infile     = q{};
my $raw_type   = q{};
my $target_dir = q{};
my $bash;

GetOptions ( 'infile=s'      => \$infile,
             'raw_type=s',   => \$raw_type,
             'target_dir=s', => \$target_dir, 
             'bash'          => \$bash, );

if ( (! $infile) or (! $raw_type) or (! $target_dir) ) { 
    die "Format: make_retro_reads_14oct2012b_rev14jan2014.pl\n",
        "    --infile|-i    [an infile listing the desired libraries]\n",
        "    --raw_types|-r [e.g., 'RNAseq' or 'jumping' ]\n",
        "    --target|-t    [target directory]\n",
        "    --bash|-b      [optionally, prepend bash header]\n",
        ;
}

if (! -e $target_dir) {
        die "Can't find required directory $target_dir\n";
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
        $desc =~ s/[.]+/./g;
        $desc =~ s/[+]/plus/g;
        $desc =~ s/\.with\.Phusion//;

        print "    mkdir $target_dir/$project ;\n";
        print "    cd $target_dir/$project ;\n";

        print '    wget --user=gec --password=gecilluminadata',
              ' --output-file ../misc/', $project, '_logfile_21may2012.txt',
              ' --recursive --level=1 --no-parent --no-directories --no-check-certificate --accept .fastq.gz',
              ' https://jumpgate.caltech.edu/runfolders/Volvox/120921_SN787_0130_BD1FPHACXX/Unaligned',
              '/Project_', $project, '_index', $index, '/Sample_', "$project ;\n",
              ;

        # For single-end data, no jumbling.  Yay!

        print '    zcat ', $project, '*_R1_*.fastq.gz > ', $target_dir, '/raw_', $raw_type, '_reads/', $project, '_', $desc, "_orig_2014.01.14.R1.fq ;\n";

        print "    cd $target_dir ;\n";

        print "    rm -rf $target_dir/$project ;\n";

        print '    gzip -9 raw_', $raw_type, '_reads/', $project, '_', $desc, "_orig_2014.01.14.R1.fq ;\n";

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

