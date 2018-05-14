#!/usr/bin/env perl

# retroname_sra_fastq_reads.pl -- Erich Schwarz <ems394@cornell.edu>, 6/20/2016; revised 8/5/2016 to add numbers to third-line headers too.
# Purpose: rename FASTQ reads from NCBI/SRS '@READ_NAME' to older Illumina '@READ_NAME#0/1[or]2', in order to unbreak ERANGE 3.1.0 [and later?] RNA-seq scripts; does keep all of header text after adding '#0/1' or '#0/2'.

# Sample input from Caltech:
# @HWI-ST0787:112:C0GPUACXX:5:1101:12825:2191 1:N:0:GATCAG
# Sample input from Berkeley:
# @HS1:159:C0JYJACXX:2:1101:2242:2150 1:N:0:
# Sample input from Edinburgh:
# @DHKW5DQ1:285:D1T8EACXX:7:1101:1397:2177 1:N:0:TATGTGGC

use strict;
use warnings;
use autodie;

my $i = 0;
my $j = 0;

my $stem        = q{};
my $distal_text = q{};
my $output      = q{};

my $input_file = $ARGV[0];
my $suffix     = $ARGV[1];

if ( ( $suffix != 1 ) and ( $suffix != 2 ) ) {
    die "Format: retroname_sra_fastq_reads.pl [input file] [suffix, either '1' or '2'] > [retro-suffixed file in STDOUT]\n";
}

open my $INPUT, '<', $input_file;
while (my $input = <$INPUT>) { 
    chomp $input;
    $i++;
    $j = ($i % 4);
    if ( ( $j != 1 ) and ( $j != 3 ) ) { 
        print "$input\n";
    }
    elsif ( $j == 1 ) { 
        if ( $input !~ /\A [@] \S+ \s .+ \z/xms ) { 
            die "Can't parse FASTQ header: $input\n";
        } 
        if ( $input =~ /\A [@] (\S+) (\s .+) \z/xms ) {
            $stem        = $1;
            $distal_text = $2;
            $output = '@' . $stem . '#0/'. $suffix . $distal_text; 
            print "$output\n";
            $stem   = q{};
            $distal_text = q{};
            $output = q{};
        }
        else {
            die "Can't parse FASTQ header: $input\n";
        }
    }
    elsif ( $j == 3 ) {
        if ( $input !~ /\A [+] \S+ \s .+ \z/xms ) {
            die "Can't parse FASTQ second-header: $input\n";
        }
        if ( $input =~ /\A [+] (\S+) (\s .+) \z/xms ) {
            $stem        = $1;
            $distal_text = $2;
            $output = '+' . $stem . '#0/'. $suffix . $distal_text;
            print "$output\n";
            $stem   = q{};
            $distal_text = q{};
            $output = q{};
        }
        else {
            die "Can't parse FASTQ header: $input\n";
        }
    }
}

