#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my @bam_files = ();

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A \S+ \. bam \z/xms ) {
        push @bam_files, $input;
    }
    else {
        die "Can't parse input line: $input\n";
    }
}

if (@bam_files) {
    print '#!/bin/bash -login', "\n";
    print '#PBS -l walltime=024:00:00', "\n";
    print '#PBS -l nodes=1:ppn=8', "\n";
    print '#PBS -l mem=32gb', "\n";
    print '#PBS -N job_pilon_nigoni_genome_2015.11.16.01.sh', "\n";
    print '#PBS -q main', "\n";
    print '#PBS -M ems394@cornell.edu', "\n";
    print '#PBS -m abe', "\n";
    print '#PBS -A ged', "\n";
    print '#PBS -r n', "\n";
    print '#PBS -V', "\n";
    print 'cd /mnt/ls15/scratch/users/emsch/work_rsync/2015/caenogens/nigoni/pilon ;', "\n";
    print 'java -Xmx31G -jar /mnt/home/emsch/src/pilon-1.14.jar --genome nigoni_mhap_pbjelly_quiver.decont_2015.11.11.fasta';
    foreach my $bam_file (@bam_files) {
        print " --frags $bam_file";
    }
    print ' --output nigoni_mhap_pbjelly_quiver.decont_2015.11.11';
    print ' --changes --chunksize 8000000 --diploid --threads 8 --verbose';
    print ' 1>nigoni_mhap_pbjelly_quiver.decont_2015.11.11.pilon_report.out.txt';
    print ' 2>nigoni_mhap_pbjelly_quiver.decont_2015.11.11.pilon_report.err.txt ;';
    print "\n";
}



