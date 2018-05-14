#!/usr/bin/env perl

use strict;
use warnings;

my @types = qw ( w_seg no_seg );

my %type2arg = ( 'w_seg'  => 'yes',
                 'no_seg' => 'no', );

my @sets = (1..10);
@sets    = map { sprintf "%02d", $_ } @sets;

foreach my $type (@types) {
    my $arg = $type2arg{$type};
    foreach my $set (@sets) {
        my $output_script = "job_psiblast_nonTreeFam_wormprots." . "$set" . '_10x_1e-04_' . "$type" . '.nematodes_and_target_TFams9.sh';
        $output_script    = safename($output_script);

        my $input_stem = 'nonTreeFam_wormprots.' . $set;
        my $input_seq  = $input_stem . '.fa';

        open my $OUTPUT, '>', $output_script or die "Can't open output script $output_script: $!";

        print $OUTPUT '#!/bin/bash -login', "\n";
        print $OUTPUT '#PBS -l walltime=110:00:00', "\n";
        print $OUTPUT '#PBS -l nodes=1:ppn=8', "\n";
        print $OUTPUT '#PBS -l mem=20gb', "\n";
        print $OUTPUT '#PBS -N ', "$output_script\n";
        print $OUTPUT '#PBS -q main', "\n";
        print $OUTPUT '#PBS -M ems394@cornell.edu', "\n";
        print $OUTPUT '#PBS -m abe', "\n";
        print $OUTPUT '#PBS -A ged-intel11', "\n";
        print $OUTPUT '#PBS -r n', "\n";
        print $OUTPUT '#PBS -V', "\n";
        print $OUTPUT "module load BLAST+/2.2.27 ;\n"; 
        print $OUTPUT "cd /mnt/home/emsch/work/DUF_R01/treefam ;\n";
        print $OUTPUT "psiblast -query $input_seq -db dbs/nematodes_and_target_TFams9",
                      " -out $input_stem", '.psiblast_10x_1e-04_', "$type", '.nematodes_and_target_TFams9.tbl.txt',
                      " -evalue 1e-04 -inclusion_ethresh 1e-04 -outfmt 6 -seg $arg -num_threads 8 -num_iterations 10 ;\n",
                      ;

        close $OUTPUT or die "Can't close filehandle to output script $output_script: $!";
    }
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

