#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;
use List::MoreUtils qw(uniq);

my $curr_dir = getcwd;

my $adaptor_file = '/mnt/ls15/scratch/users/emsch/Acey_2016/immuno_rnaseq/indiv_seqs/illumina_adaptors_all_2016.10.08.01.all.fa';

my $i = 0;

my @input_files = ();
my %file2script = ();

while (my $input = <>) {
    chomp $input;
    if ( $input !~ /\A \S+ (?: \.fastq | \.fq) \z/xms ) { 
        die "Cannot accept input file as a valid FastQ file: $input\n";
    }
    push @input_files, $input;
}

@input_files = sort @input_files;
@input_files = uniq @input_files;

foreach my $input_file (@input_files) {
    $i++;
    my $j = sprintf "%02u", $i;

    my $output_script = "job_acey_trimmomatic_count_2017.03.15.$j.sh";
    $output_script    = safename($output_script);

    $file2script{$input_file} = $output_script;
}

$i = 0;    

foreach my $input_file (@input_files) {
    my $output_file = $input_file ;

    # delete previous name-tags for filtering
    $output_file =~ s/filt1\.//g;
    $output_file =~ s/no_cont\.min40\.//g;

    # tag as Trimmo-filtered:
    $output_file =~ s/\.(fq|fastq)\z/.trim_filt.fq/;

    $output_file = safename($output_file);

    my $output_file_count = "$output_file.count.txt";
    $output_file_count    = safename($output_file_count);

    my $k = ($i + 1);

    my $output_script = $file2script{$input_file};

    open my $OUT, '>', $output_script;

    print $OUT '#!/bin/bash -login', "\n";
    print $OUT '#PBS -l walltime=003:59:00', "\n";
    print $OUT '#PBS -l nodes=1:ppn=8', "\n";
    print $OUT '#PBS -l mem=32gb', "\n";
    print $OUT "#PBS -N $output_script\n";
    print $OUT '#PBS -q main', "\n";
    print $OUT '#PBS -M ems394@cornell.edu', "\n";
    print $OUT '#PBS -m abe', "\n";
    print $OUT '#PBS -A ged', "\n";
    print $OUT '#PBS -r n', "\n";
    print $OUT '#PBS -V', "\n";
    print $OUT "cd $curr_dir ;\n";
    print $OUT 'export TRIM="/mnt/home/emsch/src/Trimmomatic-0.36" ;', "\n";
    print $OUT 'java -jar $TRIM/trimmomatic-0.36.jar SE -threads 7 -phred33 ';
    print $OUT "$input_file $output_file ";
    print $OUT 'ILLUMINACLIP:', $adaptor_file, ':2:30:10 ';
    print $OUT "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 ;\n";
    print $OUT "cat $output_file | count_simple_fastq_residues.pl > $output_file_count ;\n";

    if ( exists $input_files[$k] ) {
        my $next_output_file = $input_files[$k];
        my $next_output_script  = $file2script{$next_output_file};
        print $OUT "qsub $next_output_script ;\n";
    }

    close $OUT;
    $i++;
}

sub safename {
    my $_filename = $_[0];
    my $_orig_filename = $_filename;
    if (-e $_orig_filename) {
        my $_suffix1 = 1;
        $_filename = $_filename . ".$_suffix1";
        while (-e $_filename) {
            $_suffix1++;
            $_filename =~ s/\.\d+\z//xms;
            $_filename = $_filename . ".$_suffix1";
        }
    }
    return $_filename;
}

