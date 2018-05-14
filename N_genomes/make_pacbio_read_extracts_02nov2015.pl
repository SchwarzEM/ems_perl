#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $i = 0;

# Store very long file names here:
my $smrtwrap = '/mnt/lustre_scratch_2012/schwarz/work/2015/caenogens/smrtanalysis/install/smrtanalysis_2.3.0.140936/smrtcmds/bin/smrtwrap';
my $bash5tools  = '/mnt/lustre_scratch_2012/schwarz/work/2015/caenogens/smrtanalysis/install/smrtanalysis_2.3.0.140936/analysis/bin/bash5tools.py';

my @inputs = ();

while (my $input = <>) {
    chomp $input;

    # Ensure that the file both exists and is readable:
    if (! -r $input ) {
        die "Can't read this putative file: $input\n";
    }

    # Enforce very stereotypical input:
    # E.g.:
    # /mnt/lustre_scratch_2012/schwarz/work/2015/caenogens/wallacei/quiver/full_archives/wallacei/*/*/Analysis_Results/*.bas.h5

    if ( $input =~ /\A \/mnt\/lustre_scratch_2012\/schwarz\/work\/2015\/caenogens
                       \/wallacei\/quiver\/full_archives\/wallacei
                       \/ [^\/\s]+ \/ [^\/\s]+  \/Analysis_Results\/ [^\/\s]+ \.bas\.h5
                    \z /xms ) {
        push @inputs, $input;
    }
    else {
        die "Can't parse input: $input\n";
    }
}

my $qsub_file = "job_pacbio_read_extracts_02nov2015.sh";
$qsub_file    = safename($qsub_file);

open my $QSUB, '>', $qsub_file;

print $QSUB '#!/bin/bash -login', "\n";
print $QSUB '#PBS -l walltime=024:00:00', "\n";
print $QSUB '#PBS -l nodes=1:ppn=1', "\n";
print $QSUB '#PBS -l mem=16gb', "\n";
print $QSUB "#PBS -N $qsub_file\n";
print $QSUB '#PBS -q main', "\n";
print $QSUB '#PBS -M ems394@cornell.edu', "\n";
print $QSUB '#PBS -m abe', "\n";
print $QSUB '#PBS -A ged', "\n";
print $QSUB '#PBS -r n', "\n";
print $QSUB '#PBS -V', "\n";
print $QSUB 'cd /mnt/lustre_scratch_2012/schwarz/work/2015/caenogens/wallacei/re_extract_pacbio_reads ;', "\n";
print $QSUB "echo `date` > $qsub_file.start_time.txt ;\n";

foreach my $infile (@inputs) {
    $i++;
    my $j = sprintf "%02i", $i;
    my $outfile_prefix = "wallacei_subreads_2015.11.02.$j";

    print $QSUB "$smrtwrap $bash5tools --verbose --readType subreads --outType fastq --minLength 50 --outFilePrefix $outfile_prefix",
                " $infile",
                ' 2>', "$qsub_file.progress.$j.txt",
                " ;\n",
                ; 
    }

close $QSUB;

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

