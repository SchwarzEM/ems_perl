#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my @paired_bam_files   = ();
my @unpaired_bam_files = ();
my $genome_sequence    = 'tropicalis_hap33x.pb.qv_2016.01.12.dna.fa';
my $output_prefix      = 'tropicalis_hap33x.pb.qv.pil_2016.01.12.dna';

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A \S* unpaired \S* \. bam \z/xms ) {
        push @unpaired_bam_files, $input;
    }
    elsif ( $input =~ /\A \S* paired \S* \. bam \z/xms ) {
        push @paired_bam_files, $input;
    }
    else {
        die "Can't parse input line: $input\n";
    }
}

@paired_bam_files   = sort @paired_bam_files;
@unpaired_bam_files = sort @unpaired_bam_files;

if ( $genome_sequence and @paired_bam_files ) {
    print '#!/bin/bash -login', "\n";
    print '#PBS -l walltime=003:59:00', "\n";
    print '#PBS -l nodes=1:ppn=8', "\n";
    print '#PBS -l mem=32gb', "\n";
    print '#PBS -N job_pilon_tropicalis.33x_genome_2016.08.02.01.sh', "\n";
    print '#PBS -q main', "\n";
    print '#PBS -M ems394@cornell.edu', "\n";
    print '#PBS -m abe', "\n";
    print '#PBS -A ged', "\n";
    print '#PBS -r n', "\n";
    print '#PBS -V', "\n";
    print 'cd /mnt/ls15/scratch/users/emsch/work_rsync/2015/caenogens/tropicalis/pilon_33x ;', "\n";
    print "java -Xmx31G -jar /mnt/home/emsch/src/pilon-1.14.jar --genome $genome_sequence";
    foreach my $paired_bam_file (@paired_bam_files) {
        print " --frags $paired_bam_file";
    }
    foreach my $unpaired_bam_file (@unpaired_bam_files) {
        print " --unpaired $unpaired_bam_file";
    }
    print " --output $output_prefix";
    print ' --changes --fix bases --chunksize 10000000 --threads 8 --verbose';
    print ' 1>';
    print "$output_prefix.pilon_report.out.txt";
    print ' 2>';
    print "$output_prefix.pilon_report.err.txt ;";
    print "\n";
}

