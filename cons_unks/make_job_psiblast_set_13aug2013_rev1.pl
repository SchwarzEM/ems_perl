#!/usr/bin/env perl

use strict;
use warnings;

my @types = qw ( w_seg no_seg );

my %type2arg = ( 'w_seg'  => 'yes',
                 'no_seg' => 'no', );

my @sets = (1..100);
@sets    = map { sprintf "%03d", $_ } @sets;

foreach my $type (@types) {
    my $arg = $type2arg{$type};
    foreach my $set (@sets) {
        my $output_stem   = "job_psiblast_nonTreeFam_wormprots." . "$set" . '_10x_1e-04_' . "$type" . '.nematodes_and_target_TFams9.sh';
        my $output_script = "/mnt/home/emsch/work/DUF_R01/treefam/subscripts/". $output_stem;

        my $input_stem = 'nonTreeFam_wormprots.' . $set;
        my $input_seq  = $input_stem . '.fa';

        my $handoff  = q{};

        my $this_index = $set;
        $this_index    =~ s/\A[0]+//;

        # If we have reached set 10, 20, 30, 40, or 50, then ( $this_index % 10 ) is 0, and we do *not* want to keep going with handoffs.
        if ( $this_index % 10 ) {
            my $next_index = $this_index;
            $next_index++;
            my $next_set = sprintf "%03d", $next_index;
            my $next_script = "/mnt/home/emsch/work/DUF_R01/treefam/subscripts/job_psiblast_nonTreeFam_wormprots." 
                                   . "$next_set" . '_10x_1e-04_' . "$type" . '.nematodes_and_target_TFams9.sh';
            $handoff  = "qsub $next_script ;";
        }

        open my $OUTPUT, '>', $output_script or die "Can't open output script $output_script: $!";

        print $OUTPUT '#!/bin/bash -login', "\n";
        print $OUTPUT '#PBS -l walltime=011:00:00', "\n";
        print $OUTPUT '#PBS -l nodes=1:ppn=8', "\n";
        print $OUTPUT '#PBS -l mem=20gb', "\n";
        print $OUTPUT '#PBS -N ', "$output_stem\n";
        print $OUTPUT '#PBS -q main', "\n";
        print $OUTPUT '#PBS -M ems394@cornell.edu', "\n";
        print $OUTPUT '#PBS -m abe', "\n";
        print $OUTPUT '#PBS -A ged-intel11', "\n";
        print $OUTPUT '#PBS -r n', "\n";
        print $OUTPUT '#PBS -V', "\n";
        print $OUTPUT "module load BLAST+/2.2.27 ;\n"; 
        print $OUTPUT "cd /mnt/home/emsch/work/DUF_R01/treefam ;\n";
        print $OUTPUT "psiblast -query subseqs/", "$input_seq -db dbs/nematodes_and_target_TFams9",
                      " -out subblasts/", "$input_stem", '.psiblast_10x_1e-04_', "$type", '.nematodes_and_target_TFams9.tbl.txt',
                      " -evalue 1e-04 -inclusion_ethresh 1e-04 -outfmt 6 -seg $arg -num_threads 8 -num_iterations 10 ;\n",
                      ;

        # We only append this *if* it has been defined, i.e., if we haven't reached 10 etc.
        print $OUTPUT "$handoff\n" if $handoff;

        close $OUTPUT or die "Can't close filehandle to output script $output_script: $!";
    }
}

sub increment_serial_no {
    my $_serial_no = $_[0];
    $_serial_no =~ s/\A[0]+//;
    my $_next_no = $_serial_no;
    $_next_no++;
    return $_next_no++;
}

