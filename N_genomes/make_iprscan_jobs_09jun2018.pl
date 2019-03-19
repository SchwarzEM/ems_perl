#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Spec::Functions;  # catdir
use Scalar::Util qw(looks_like_number);

my $stem    = q{};
my $workdir = q{};
my $limit   = q{};
my $censor  = q{};

$stem    = $ARGV[0] if $ARGV[0];
$workdir = $ARGV[1] if $ARGV[1];
$limit   = $ARGV[2] if $ARGV[2];
$censor  = $ARGV[3] if $ARGV[3];

if ( ( $stem !~ /\A \S+ \z/xms ) or (! looks_like_number($limit) ) or ( $limit < 1 ) or ( $limit != int($limit) ) ) {
    die "make_iprscan_jobs_01dec2016.pl\n",
        "    [stem -- must be single block of text]\n",
        "    [planned working directory for the scripts]\n",
        "    [limit -- must be positive integer]\n",
        "    [optional censored qsub jobs, listed as X,Y,Z -- joined with commas]\n",
        ;
}

my @no_qsub_nos = ();
my %no_qsub     = ();

if ($censor) {
    @no_qsub_nos = split q{,}, $censor;
    foreach my $banned (@no_qsub_nos) {
        $banned = sprintf "%02i", $banned;
        $no_qsub{$banned} = 1;
    }
}

foreach my $i (1..$limit) {
    my $file_no      = sprintf "%02i", $i;
    my $next_file_no = sprintf "%02i", ($i + 1);

    my $input_file = "$stem.$file_no.fa";
    if (! -r $input_file ) {
        die "Cannot read putative input file: $input_file\n";
    }
    # Convert this to a full pathname for iprscan:
    $input_file = catfile($workdir, $input_file);

    my $output_stem = "$stem.$file_no";
    # again, make into a full path:
    $output_stem = catfile($workdir, $output_stem);

    my $job1 = "job_iprscan_2018.06.09.$file_no.sh";
    my $job2 = "job_iprscan_2018.06.09.$next_file_no.sh";

    $job1 = safename($job1);
    $job2 = safename($job2);

    open my $JOB1, '>', $job1 ;

    print $JOB1 '#!/bin/bash -login', "\n";
    print $JOB1 '#PBS -l walltime=003:59:00', "\n";
    print $JOB1 '#PBS -l nodes=1:ppn=5', "\n";
    print $JOB1 '#PBS -l mem=32gb', "\n";
    print $JOB1 "#PBS -N $job1\n";
    print $JOB1 '#PBS -q main', "\n";
    print $JOB1 '#PBS -M ems394@cornell.edu', "\n";
    print $JOB1 '#PBS -m abe', "\n";
    print $JOB1 '#PBS -A ged', "\n";
    print $JOB1 '#PBS -r n', "\n";
    print $JOB1 '#PBS -V', "\n";
    print $JOB1 "cd $workdir ;\n";
    print $JOB1 "module load InterProScan/5.18-57.0 ;\n";
    print $JOB1 "interproscan.sh -dp -hm -iprlookup -goterms -i $input_file -b $output_stem -T temp.$file_no ;\n";

    if (! exists $no_qsub{$next_file_no} ) { 
        print $JOB1 "qsub $job2 ;\n";
    }

    close $JOB1;
}

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

