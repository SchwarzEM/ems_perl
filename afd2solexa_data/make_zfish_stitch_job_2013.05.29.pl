#!/usr/bin/env perl

use strict;
use warnings;

my $temp         = 'temp_work_file_2013.05.29.txt';
my $final_output = 'zfish_reads_file_2013.05.29.pe.fq';
$final_output    = safename($final_output);

my @work_steps   = ();

while (my $file1 = <>) { 
    chomp $file1;
    if ( $file1 =~ /\A (\S+_)1(\.fastq\.gz) \z/xms ) { 
        my $text1 = $1;
        my $text2 = $2;
        my $file2 = $text1 . '2' . $text2 ;
        if ( (! -e $file1 ) or (! -e $file2 ) ) { 
            die "Can't find both of these files: $file1; $file2\n";
        }
        my $work_step = q{};
        $work_step = "rm $temp ;\n";
        push @work_steps, $work_step;
        $work_step = q{shuffleSequences_fastq.pl <(zcat } . $file1 . q{) <(zcat } . $file2 . ") $temp ;\n";
        push @work_steps, $work_step;
        $work_step = "cat $temp >> $final_output ;\n";
        push @work_steps, $work_step;
    }
}

print '#!/bin/bash -login', "\n";
print '#PBS -l walltime=024:00:00', "\n";
print '#PBS -l nodes=1:ppn=1', "\n";
print '#PBS -l mem=10gb', "\n";
print '#PBS -N job_stitch_all_zfish_rseq_data_2013.05.29', "\n";
print '#PBS -q main', "\n";
print '#PBS -M ems394@cornell.edu', "\n";
print '#PBS -m abe', "\n";
print '#PBS -A ged-intel11', "\n";
print '#PBS -r n', "\n";
print '#PBS -V', "\n";
print 'cd /mnt/scratch/emsch/ALA ;', "\n";
print "touch $temp ;\n";
print @work_steps;
print "rm $temp ;\n";
print q{/mnt/home/emsch/perl.svn/trunk/illumina/count_simple_fastq_residues.pl }, $final_output, q{ > }, $final_output, ".count.txt ;\n";

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

