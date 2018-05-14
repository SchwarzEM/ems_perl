#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;

my $working_dir = getcwd;

my %seen = ();

while (my $input = <>) {
    chomp $input;
    if ( (! -r $input ) or ( $input !~ /\A \S+ \z/xms ) ) {
        die "Can't read this putative file: $input\n";
    }
    $seen{$input} = 1;
}

my @inputs = sort keys %seen;

print '#!/bin/bash -login', "\n";
print '#PBS -l walltime=003:59:00', "\n";
print '#PBS -l nodes=1:ppn=8', "\n";
print '#PBS -l mem=64gb', "\n";
print "#PBS -N job_orphan_blastp_DATE.sh\n";
print '#PBS -q main', "\n";
print '#PBS -M ems394@cornell.edu', "\n";
print '#PBS -m abe', "\n";
print '#PBS -A ged', "\n";
print '#PBS -r n', "\n";
print '#PBS -V', "\n";
print "cd $working_dir ;\n";
print "module load BLAST+/2.2.31 ;\n";

foreach my $input (@inputs) {
    my $output = "$input.blastp_E.1e-03.no_seg.orthfind_5spp_2016.11.05.proteome.tsv.txt";
    $output    = safename($output);

    print "blastp -num_threads 8 -max_target_seqs 1000000 -outfmt '6 qseqid sseqid evalue'";
    print " -evalue 1e-03 -seg no ";
    print " -db /mnt/ls15/scratch/users/emsch/2016/caenogens/orthofinder/dbs/orthfind_5spp_2016.11.05.01.proteome";
    print " -query $input -out $output ;\n";
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

