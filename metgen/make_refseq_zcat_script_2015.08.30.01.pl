#!/usr/bin/env perl

# make_refseq_zcat_script_2015.08.30.01.pl -- Erich Schwarz <ems394@cornell.edu>, 8/30/2015.
# Purpose: given many 'nonredundant' RefSeq *.faa.gz files with non-Unixy numbering, make it easy to make a shell script that zcats them in order: 1, 2, 3, ... 9, 10, ... etc.

use strict;
use warnings;
use autodie;

my $target_file = q{};
my $file_list   = `ls * | more`;
my @orig_files  = split "\n", $file_list;
my %file2num    = ();

$target_file = $ARGV[0] if $ARGV[0];

if ( $target_file !~ /\A \S+ \z/xms ) {
    die "Must provide a target file name as argument: make_refseq_zcat_script_2015.08.30.01.pl [all_non_space_arg]\n";
}

foreach my $file (@orig_files) {
    if ( $file =~ /\A \S+ \.nonredundant_protein\. (\d+) \.protein\.faa\.gz \z/xms ) {
        my $num = $1;
        $file2num{$file} = $num;
    }
}

my @filt_files = sort { $file2num{$a} <=> $file2num{$b} } keys %file2num;

my $header      = '#!/bin/bash' . "\n\n";
my $text_action = '> ';

foreach my $file (@filt_files) {
    # Print the header, once:
    print $header if $header;
    $header = q{};

    # For the very first line, use '> ' so that older versions of $target_file are erased:
    print "    zcat $file $text_action $target_file ;\n";

    # For every succeeding line, catenate more text nondestructively with '>>':
    $text_action = '>>';
}

