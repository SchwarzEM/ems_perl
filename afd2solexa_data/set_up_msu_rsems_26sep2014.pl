#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use File::Basename;
use File::Spec;

my $read_table = q{};
my $work_dir   = q{};
my $tag        = q{};

$read_table = $ARGV[0] if $ARGV[0];
$work_dir   = $ARGV[1] if $ARGV[1];
$tag        = $ARGV[2] if $ARGV[2];   # e.g. '25sep2014'

my %qual_arguments = (
    'Sanger FASTQ format'                           => '--phred33-quals',  # Input quality scores are encoded as Phred+33. (Default in RSEM: on)
    'Illumina FASTQ format, Illumina pipeline 1.3+' => '--phred64-quals',  # Input quality scores are encoded as Phred+64 (default for GA Pipeline ver. >= 1.3). (Default: off)
);

# Always require that a new directory be created for jobs.
if ( (-e $work_dir) or ( $work_dir !~ /\A \S+ \z/xms ) ) {
    warn "Usage: set_up_msu_rsems_26sep2014.pl [read table, with infiles and qual types] [new working directory] [tag for scripts/subdirectories]\n";
    die "Preexisting work directory: \"$work_dir\"\n";
}
if ( $work_dir =~ /\A \S+ \z/xms ) {
    mkdir $work_dir, 0755;
}

if (! $tag) { 
    warn "Usage: set_up_msu_rsems_24sep2014.pl [read table, with infiles and qual types] [new working directory] [tag for scripts/subdirectories]\n";
    die "Please specify some tag for scripts and subdirectories, e.g., \"25sep2014\"\n";
}

open my $TABLE, '<', $read_table;
while (my $input = <$TABLE>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t ([^\t]+) \z/xms ) {
        my $infile = $1;
        my $qual   = $2;

        if (! exists $qual_arguments{$qual} ) {
            warn "Usage: set_up_msu_rsems_24sep2014.pl [read table, with infiles and qual types] [working directory] [tag for scripts/subdirectories]\n";
            die "From read table $read_table, can't parse quality type in input line: $input\n";
        }
        if (! -e $infile) {
            warn "Usage: set_up_msu_rsems_24sep2014.pl [read table, with infiles and qual types] [working directory] [tag for scripts/subdirectories]\n";
            die "From read table $read_table, nonexistent input file in input line: $input\n";
        }

        my $qual_argument = $qual_arguments{$qual};

        my $basename = basename $infile;
        $basename =~ s/\.fq\z//;
        $basename =~ s/\.fastq\z//;
        $basename =~ s/\.fa\z//;
        $basename =~ s/\.fasta\z//;

        my $work_subdir = $basename . '_' . $tag . '_dir';
        $work_subdir = File::Spec->catfile($work_dir, $work_subdir);
        $work_subdir = safename($work_subdir);

        mkdir $work_subdir, 0755;

        my $work_label      = 'rsem_' . $basename . '_' . $tag;
        my $job_script_name = 'job_' . $work_label;

        my $job_script      = $job_script_name;
        $job_script         = $job_script . '.sh';
        $job_script         = File::Spec->catfile($work_subdir, $job_script);
        $job_script         = safename($job_script);

        open my $JOB_SCRIPT, '>', $job_script;

        print $JOB_SCRIPT '#!/bin/bash -login', "\n";
        print $JOB_SCRIPT '#PBS -l walltime=008:00:00', "\n";
        print $JOB_SCRIPT '#PBS -l nodes=1:ppn=8', "\n";
        print $JOB_SCRIPT '#PBS -l mem=21gb', "\n";
        print $JOB_SCRIPT "#PBS -N $job_script_name\n";
        print $JOB_SCRIPT '#PBS -q main', "\n";
        print $JOB_SCRIPT '#PBS -M ems394@cornell.edu', "\n";
        print $JOB_SCRIPT '#PBS -m abe', "\n";
        print $JOB_SCRIPT '#PBS -A ged', "\n";
        print $JOB_SCRIPT '#PBS -r n', "\n";
        print $JOB_SCRIPT '#PBS -V', "\n";
        print $JOB_SCRIPT "cd $work_subdir ;\n";
        print $JOB_SCRIPT "rsem-calculate-expression --bowtie2 -p 8 --no-bam-output --calc-pme $qual_argument";
        print $JOB_SCRIPT ' --fragment-length-mean 200 --fragment-length-sd 20 --estimate-rspd --ci-memory 20000';
        print $JOB_SCRIPT " $infile /mnt/home/emsch/work/AFD_etc/rsem_data/c_elegans.WS245_all.txs.rnaseq $work_label ;\n";

        close $JOB_SCRIPT;
    }
    else {
        die "From read table $read_table, can't parse input line: $input\n";
    }
}
close $TABLE;

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


