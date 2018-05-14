#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use File::Basename;

my @input_files   = @ARGV;
my $temp_read_set = q{};
my $header        = '#!/bin/bash' . "\n\n";

foreach my $input_file (@input_files) {
    my $base_input_file        = basename $input_file;
    my $temp_sample_reads_file = "$base_input_file.tmp";
    $temp_sample_reads_file    = safename($temp_sample_reads_file);

    print $header if $header;
    $header = q{};

    print "    echo \"Analyzing RNA-seq file: $input_file\" ;\n";
    print "    head --lines=40000 $input_file > $temp_sample_reads_file ;\n";
    print "    SolexaQA++ analysis $temp_sample_reads_file ;\n";
    print "    rm $temp_sample_reads_file.* ;\n";
    print "    echo ;\n";
    print "\n";
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

