#!/usr/bin/env perl

# retrofix_fastq_reads.pl -- Erich Schwarz <ems394@cornell.edu>, 8/7/2021.
# Purpose: fix FASTQ reads from quasi-half-retro "/1", "/2" tags to fully retro Illumina '@READ_NAME#0/1'"

use strict;
use warnings;
use autodie;

use Getopt::Long;

# Sample input from ENA/NCBI:
#
# @SRR2125612.2001 HWI-ST514:143982632:C37PRACXX:5:1101:19721:2677/1
# @SRR2125612.2001 HWI-ST514:143982632:C37PRACXX:5:1101:19721:2677/2
#
# @SRR2138602.6846449 6846449/1
# @SRR2138602.6846449 6846449/2 

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
            if ( $input !~ /\A [@] \S+ \s \S .+ \/ [12] \s* \z/xms ) { 
                die "Can't parse FASTQ header: $input\n";
            } 
            if ( $input =~ /\A [@] (\S+) (\s \S .+ \/ ([12]) \s*) \z/xms ) {
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

