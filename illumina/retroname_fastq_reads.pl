#!/usr/bin/env perl

# retroname_fastq_reads.pl -- Erich Schwarz <ems394@cornell.edu>, 11/22/2019.
# Purpose: rename FASTQ reads from newer Illumina '@READ_NAME 1:N:0:GATCAG' to older Illumina '@READ_NAME#0/1', in order to unbreak ERANGE 3.1.0 [and later?] RNA-seq scripts; does keep all of header text after adding '#0/1' or '#0/2'; has '--racon' option to append '1' or '2'.

use strict;
use warnings;
use autodie;

use Getopt::Long;

# Sample input from Caltech:
# @HWI-ST0787:112:C0GPUACXX:5:1101:12825:2191 1:N:0:GATCAG
# Sample input from Berkeley:
# @HS1:159:C0JYJACXX:2:1101:2242:2150 1:N:0:
# Sample input from Edinburgh:
# @DHKW5DQ1:285:D1T8EACXX:7:1101:1397:2177 1:N:0:TATGTGGC

my @infiles = ();

my $i = 0;
my $j = 0;

my $stem        = q{};
my $digit       = q{};
my $distal_text = q{};
my $output      = q{};

my $racon;
my $help;

GetOptions ( 'infiles=s{,}' => \@infiles,
             'racon'        => \$racon,
             'help'         => \$help,   );

if ( $help or (! @infiles) ) { 
    die "Format: retroname_fastq_reads.pl\n",
        "    --infile|-i  <input stream/files>\n",
        "    --racon|-r   [add Racon-friendly suffixes '1' or '2' instead of retro-suffixes]\n",
        "    --help|-h    [print this message]\n",
        ;
}

foreach my $infile (@infiles) {
    my $INPUT_FILE;
    if ($infile eq '-') {
        # Special case: get the stdin handle
        $INPUT_FILE = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $INPUT_FILE, '<', $infile or die "Can't open input file $infile. $!\n";
    }
    while (my $input = <$INPUT_FILE>) { 

        chomp $input;
        $i++;
        $j = ($i % 4);
        if ( $j != 1 ) { 
            print "$input\n";
        }
        if ( $j == 1 ) { 
            if ( $input !~ /\A [@] \S+ \s \d:(?:Y|N):\d+:[ACGTN+]* \s* \z/xms ) { 
                die "Can't parse FASTQ header: $input\n";
            } 
            if ( $input =~ /\A [@] (\S+) (\s (\d):(?:Y|N):\d+:[ACGTN+]* \s*) \z/xms ) {
                $stem        = $1;
                $distal_text = $2;
                $digit       = $3;
                my $suffix   = q{};

                if ($racon) {
                    $suffix = $digit;
                }
                else {
                    $suffix = '#0/'. $digit;
                }
                $output = '@' . $stem . $suffix . $distal_text; 
                print "$output\n";
                $stem   = q{};
                $digit  = q{};
                $output = q{};
            }
        }
    }
}

