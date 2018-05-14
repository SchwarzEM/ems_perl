#!/usr/bin/env perl

# make_StatGap_tsv.pl -- Erich Schwarz <emsch@@its.caltech.edu>, 3/2/2013.
# Purpose: given a TSV or a slice of a TSV, make a table which can be either the 'vec' or the 'sim' input for Statistics::Gap.

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my @input_files = ();
my $sim;
my $vec;
my $help;

my $tab_count   = 0;
my @input_lines = ();

GetOptions ( 'input_files=s{,}' => \@input_files,
             'sim'              => \$sim,
             'vec'              => \$vec,
             'help'             => \$help, );

if ( $help or (! @input_files ) or ( $sim and $vec ) or ( (! $sim) and (! $vec) ) ) {
    die "Format: make_StatGap_tsv.pl\n",
        "        --input_files|-i        [input file(s), or '-' if stream]\n",
        "        --sim|-s [or] --vec|-v  [mutually exclusive; for most purposes, will want 'vec'; either way, dense matrix format]\n",
        "        --help|-h               [print this message]\n",
        ;
}

# Accept either a stream from '-' or a standard file.
my $INPUT_FILE;
foreach my $infile (@input_files) { 
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

        # Enforce consistent tab counts for all lines.
        my $line_tab_count = ( $input =~ tr/\t/\t/ );
        if ( $line_tab_count != $tab_count ) {
            if ( $tab_count > 0 ) {
                die "Inconsistent tab counts (first $tab_count, then $line_tab_count)\n";
            }
            $tab_count = $line_tab_count;
        }

        # Ensure that every single value is, in fact, a number:
        my @vals = split "\t", $input;
        foreach my $val (@vals) { 
            if (! looks_like_number($val)) { 
                die "Non-numerical value $val in input line: $input\n";
            }
        }

        # Store line, if it is really OK.
        push @input_lines, $input;
    }
    close $INPUT_FILE or die "Can't close filehandle to input file $infile: $!\n";
}

my $line_count = @input_lines;

# Remember, there are +1 more columns than tabs.
$tab_count++;

# Now that both values exist and are corrected, enforce their equality for 'sim'-compatible matrix outputs.
if ( $sim and ( $line_count != $tab_count ) ) {
    die "Cannot have 'sim'-type matrix output from input whose line count ($line_count) and column count ($tab_count) are not equal.\n";
}

# For 'vec' input, Statistics::Gap wants #rows #cols header line.
# For 'sim' input, Statistics::Gap wants just #rows header line.

print "$line_count $tab_count\n" if $vec;
print "$line_count\n" if $sim;

foreach my $input_line (@input_lines) { 
    print "$input_line\n";
}

