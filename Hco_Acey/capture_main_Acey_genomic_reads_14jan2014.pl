#!/usr/bin/env perl

# capture_main_Acey_genomic_reads_14jan2014.pl -- Erich Schwarz <ems394@cornell.edu>, 1/14/2014.

use strict;
use warnings;
use Getopt::Long;

my $infile     = q{};
my $source     = q{};
my $raw_type   = q{};
my $target_dir = q{};
my $bash;

GetOptions ( 'infile=s'      => \$infile,
             'source=s'      => \$source,
             'raw_type=s',   => \$raw_type,
             'target_dir=s', => \$target_dir, 
             'bash'          => \$bash, );

if ( (! $infile) or (! $source) or (! $raw_type) or (! $target_dir) ) { 
    die "Format: capture_main_Acey_genomic_reads_14jan2014.pl\n",
        "    --infile|-i    [an infile listing the desired libraries]\n",
        "    --source|-i    [source directory of the desired libraries]\n",
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
    # Lane #7 : (12443) Index #11 ACEY HWORM 15SEP11 (None)

    if ( $input =~ /\A Lane [ ] [#] \d+ [ ] : [ ] \( (\d+) \) [ ] Index [ ] [#] (\d+) [ ] (\S.+\S) [ ] \( None \) \s* \z/xms ) { 
        my $project = $1;
        my $index   = $2;
        my $desc    = $3;

        $desc =~ s/\[A\.[ ]ceylanicum\]/ACEY/g;
        $desc =~ s/\s/./g;
        $desc =~ s/\.with\.Phusion//;
        $desc =~ s/.15SEP11//;
        $desc =~ s/HWORM/HWORM_GENOMIC/;

        print "    mkdir $target_dir/$project ;\n";
        print "    cd $target_dir/$project ;\n";

        print '    zcat ', "$source/$project", '*_R1_*.fastq.gz > ', $target_dir, '/raw_', $raw_type, '_reads/', $project, '_', $desc, "_orig_2014.01.14.R1.fq ;\n";

        print '    zcat ', "$source/$project", '*_R2_*.fastq.gz > ', $target_dir, '/raw_', $raw_type, '_reads/', $project, '_', $desc, "_orig_2014.01.14.R2.fq ;\n";

        print "    cd $target_dir ;\n";

        print "    rm -rf  $target_dir/$project ;\n";

        print '    shuffleSequences_fastq.pl',
              ' raw_', $raw_type, '_reads/', $project, '_', $desc, '_orig_2014.01.14.R1.fq',
              ' raw_', $raw_type, '_reads/', $project, '_', $desc, '_orig_2014.01.14.R2.fq',
              ' raw_', $raw_type, '_reads/', $project, '_', $desc, "_orig_2014.01.14.fq ;\n",
              ; 

        print '    rm raw_', $raw_type, '_reads/', $project, '_', $desc, '_orig_2014.01.14.R1.fq raw_', $raw_type, '_reads/', $project, '_', $desc, '_orig_2014.01.14.R2.fq', " ;\n";

        print '    gzip -9 raw_', $raw_type, '_reads/', $project, '_', $desc, "_orig_2014.01.14.fq ;\n";

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

