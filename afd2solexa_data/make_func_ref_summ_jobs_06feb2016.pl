#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $read_file  = 0;
my $error_file = safename('func_ref_summ_jobs_06feb2016.err');
my $header     = '#!/bin/bash';

while (my $input = <>) {
    if ( $input =~ /See [ ] results [ ] in [ ] files:\s* \z/xms ) { 
        $read_file = 1;
    }
    elsif ( $read_file and ( $input =~ /\A [ ][-][ ] (\S+) /xms ) ) { 
        my $input_file = $1;
        $input_file =~ s/\\//g;
        my $output_file = q{};

        # Sample input file:
        # /mnt/home/emsch/work/2015/adrienne/go_analysis/atml1-3.vs.Col_WT.downreg.2016.02.06.func_bool.analysis/refinement-*-0.01_0.01.txt
        if ( $input_file =~ /\A (\S+) \/refinement-[*]-0.01_0.01.txt \z/xms ) { 
            my $stem = $1;
            $output_file = "$stem.tsv.txt";
            $output_file = safename($output_file);
        }
        else {
            die "Cannot parse input file: $input_file\n";
        }

        $read_file = 0;

        print "$header\n\n" if $header;
        $header = q{};

        print "    filter_func_refinement_01jun2015.pl -o -p 0.01 -i $input_file 1>";
        print "$output_file 2>>";
        print "$error_file ;\n";
    }
}

print "\n" if (! $header);

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

