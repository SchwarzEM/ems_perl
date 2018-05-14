#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my @queries = <>;
chomp @queries;

my %query2output = ();

foreach my $query (@queries) {
    if ( $query =~ /\A \.\.\/ indiv_seqs \/ (\S+) \.fa \z/xms ) { 
        my $seqname = $1;
        my $output = "$seqname.blastp_1e-06.titus.prairie_cdhit0.90.txt";
        $query2output{$query} = $output;
    }
    else {
        die "Can't parse query: $query\n";
    }
}

print '#!/bin/bash -login', "\n";
print '#PBS -l walltime=024:00:00', "\n";
print '#PBS -l nodes=1:ppn=08', "\n";
print '#PBS -l mem=16gb', "\n";
print '#PBS -N job_blastp_titus.prairie_2015.05.10.01', "\n";
print '#PBS -q main', "\n";
print '#PBS -M ems394@cornell.edu', "\n";
print '#PBS -m abe', "\n";
print '#PBS -A ged', "\n";
print '#PBS -r n', "\n";
print '#PBS -V', "\n";
print 'cd /mnt/lustre_scratch_2012/schwarz/work/2015/fourth_domain/soil_mgen ;', "\n";
print 'module load BLAST+/2.2.30 ;', "\n";

foreach my $query (@queries) {
    print "blastp -num_threads 8 -seg yes",
          " -db dbs/iowa_prairie.megahit_0.2.0_10may2015",
          " -evalue 1e-06 -num_descriptions 50 -num_alignments 50 -query $query -out $query2output{$query} ;",
          "\n",
          ;

}


