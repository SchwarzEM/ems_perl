#!/usr/bin/env perl

# light_read_sizefilter.pl -- Erich Schwarz <emsch@its.caltech.edu>, 5/31/2011.
# Purpose: given a huge set of Illumina reads in *single-line* FASTA format, do fast size-filtering without heavy load on memory.

use strict;
use warnings;
use Getopt::Long;

my @input_files = ();
my $header      = q{};
my $sequence    = q{};
my $seq_length  = 0;
my $min_length  = 0;

my $reading_header = 1;
my $help;

GetOptions ( 'input_files=s{,}' => \@input_files,
             'min=s'            => \$min_length,
             'help'             => \$help, );

if ( $help or (! @input_files ) or (! $min_length ) ) { 
    die "\n",
        "Format: light_read_sizefilter.pl\n",
        "        --input_files|-i  [input file(s), or '-' if stream]\n",
        "        --min|-m          [minimum allowable read length]\n",
        "        --help|-h\n",
        "\n",
        ;
}

# Accept either a stream from '-' or a standard file.
foreach my $infile (@input_files) { 
    my $INPUT_FILE;
    if ($infile eq '-') {
        # Special case: get the stdin handle
        $INPUT_FILE = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $INPUT_FILE, '<', $infile or die "Can't open input file $infile. $!\n";
    }

    # Size-filter the incoming FASTA reads, source by source:
    while (my $input = <$INPUT_FILE>) {
        chomp $input;
        if ( $reading_header and (! $header) and ( $input =~ /\A (> \S+ .*) \z/xms ) ) { 
            $header = $1;
            $reading_header = 0;
        }
        elsif ( (! $reading_header) and $header and ( $input =~ /\A ([ACGT]+) \z/xms ) ) { 
            $sequence = $1;
            my $seq_length = length($sequence);

            # Print qualified read headers and sequence single-lines:
            if ( $seq_length >= $min_length ) { 
                print "$header\n";
                print "$sequence\n";
            }

            # Zero stored data and seqlength:
            $header     = q{};
            $sequence   = q{};
            $seq_length = 0;

            # Get ready to read another header:
            $reading_header = 1
        }
        else { 
            die "Can't parse input line: $input\n";
        }
    }
    close $INPUT_FILE or die "Can't close filehandle to input file $infile: $!\n";
}

