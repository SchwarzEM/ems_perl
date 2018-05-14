#!/usr/bin/env perl

use strict;
use warnings;

my %stem_dirs = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (\S+)\.fasta\.\d+\.dir \z/xms ) { 
        my $stem = $1;
        $stem_dirs{$stem} = 1;
    } 
    elsif ( $input =~ /\A (\S+)\.fa\.\d+\.dir \z/xms ) {
        my $stem = $1;
        $stem_dirs{$stem} = 1;
    }
    elsif ( $input =~ /\A (\S+)\.stragglers_98i\.dir \z/xms ) {
        my $stem = $1;
        $stem_dirs{$stem} = 1;
    }
    else { 
        die "Can't parse input line: $input\n";
    }
}

print '#!/bin/bash', "\n\n";

my @stems = sort keys %stem_dirs;
foreach my $stem (@stems) { 
    my $output_file = $stem . '_total_iprscan_raw.txt';
    $output_file = safename($output_file);
    print '    cut -f 4-6,9 ', $stem, '.*.dir/*.iprscan.raw', " > $output_file;\n";
}

print "\n";

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

