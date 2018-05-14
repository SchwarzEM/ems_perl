#!/usr/bin/env perl

# retroname_fastq_reads.pl -- Erich Schwarz <emsch@caltech.edu>, 7/17/2014.
# Purpose: rename FASTQ reads from newer Illumina '@READ_NAME 1:N:0:GATCAG' to older Illumina '@READ_NAME#0/1', in order to unbreak ERANGE 3.1.0 [and later?] RNA-seq scripts; does keep all of header text after adding '#0/1' or '#0/2'.

# Sample input from Caltech:
# @HWI-ST0787:112:C0GPUACXX:5:1101:12825:2191 1:N:0:GATCAG
# Sample input from Berkeley:
# @HS1:159:C0JYJACXX:2:1101:2242:2150 1:N:0:
# Sample input from Edinburgh:
# @DHKW5DQ1:285:D1T8EACXX:7:1101:1397:2177 1:N:0:TATGTGGC

use strict;
use warnings;

my $i = 0;
my $j = 0;

my $stem        = q{};
my $digit       = q{};
my $distal_text = q{};
my $output      = q{};

while (my $input = <>) { 
    chomp $input;
    $i++;
    $j = ($i % 4);
    if ( $j != 1 ) { 
        print "$input\n";
    }
    if ( $j == 1 ) { 
        if ( $input !~ /\A [@] \S+ \s \d:(?:Y|N):\d+:[ACGTN]* \s* \z/xms ) { 
            die "Can't parse FASTQ header: $input\n";
        } 
        if ( $input =~ /\A [@] (\S+) (\s (\d):(?:Y|N):\d+:[ACGTN]* \s*) \z/xms ) {
            $stem        = $1;
            $distal_text = $2;
            $digit       = $3;
            $output = '@' . $stem . '#0/'. $digit . $distal_text; 
            print "$output\n";
            $stem   = q{};
            $digit  = q{};
            $output = q{};
        }
    }
}

