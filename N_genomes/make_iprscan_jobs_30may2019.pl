#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Spec::Functions;  # catdir
use Scalar::Util qw(looks_like_number);
use Cwd;

my $stem   = q{};
my $limit  = q{};
my $censor = q{};

$stem   = $ARGV[0] if $ARGV[0];
$limit  = $ARGV[1] if $ARGV[1];
$censor = $ARGV[2] if $ARGV[2];

my $workdir = getcwd;

if ( ( $stem !~ /\A \S+ \z/xms ) or (! looks_like_number($limit) ) or ( $limit < 1 ) or ( $limit != int($limit) ) ) {
    die "make_iprscan_jobs_30may2019.pl\n",
        "    [stem -- must be single block of text]\n",
        "    [limit -- must be positive integer]\n",
        "    [optional censored sbatch jobs, listed as X,Y,Z -- joined with commas]\n",
        ;
}

my @no_sbatch_nos = ();
my %no_sbatch     = ();

if ($censor) {
    @no_sbatch_nos = split q{,}, $censor;
    foreach my $banned (@no_sbatch_nos) {
        $banned = sprintf "%02i", $banned;
        $no_sbatch{$banned} = 1;
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

    my $job1 = "job_iprscan_Acey_2019.05.30.$file_no.sh";
    my $job2 = "job_iprscan_Acey_2019.05.30.$next_file_no.sh";

    $job1 = safename($job1);
    $job2 = safename($job2);

    open my $JOB1, '>', $job1 ;

    print $JOB1 '#!/bin/bash', "\n";
    print $JOB1 '#SBATCH --nodes=1', "\n";
    print $JOB1 '#SBATCH --partition=RM-shared', "\n";
    print $JOB1 '#SBATCH --time=04:00:00', "\n";
    print $JOB1 '#SBATCH --ntasks-per-node=8', "\n";
    print $JOB1 '#SBATCH --constraint=EGRESS', "\n";
    print $JOB1 '#SBATCH --job-name=', "$job1\n";
    print $JOB1 '#SBATCH --mail-type=ALL', "\n";
    print $JOB1 "cd $workdir ;\n";
    print $JOB1 '$SCRATCH/InterPro/interproscan-5.35-74.0/interproscan.sh';
    print $JOB1 " -cpu 7 -dp -iprlookup -goterms -i $input_file -b $output_stem -T temp.$file_no ;\n";

    if (! exists $no_sbatch{$next_file_no} ) { 
        print $JOB1 "sbatch $job2 ;\n";
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

