#!/usr/bin/env perl

use strict ;
use warnings ;
use autodie ;

my $rnaseq_file  = q{};
my $quality_type = q{};
my $report_text  = q{};

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A Analyzing [ ] RNA-seq [ ] file: [ ] (\S+) \s* \z/xms ) { 
        $rnaseq_file = $1;
    }
    elsif ( $input =~ /\A Automatic [ ] format [ ] detection: [ ] (\S.+\S) \s* \z/xms ) {
        $quality_type = $1;
        my $report_line = "$rnaseq_file\t$quality_type\n";
        if ($rnaseq_file) {
            $report_text   .= $report_line;
        }
        $rnaseq_file  = q{};
        $quality_type = q{};
    }
}

print $report_text;
