#!/usr/bin/env perl

# fa_numb_names.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/7/2010; updated 8/4/2016.
# Purpose: replace FASTA names with '0000x' series; goes well with tag_FASTA_names.pl; can optionally start with number of 2+.

use strict;
use warnings;
use autodie;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my $infile = q{};
my $header = q{};
my $pad    = 0;
my $i      = 1;
my $index  = q{};
my $help;

GetOptions ( 'infile=s'  => \$infile,
             'padding=s' => \$pad,
             'start=i'   => \$i,
             'help'      => \$help
           );

if ( $help or (! $infile) ) { 
    die "Format: fa_numb_names.pl\n",
        "    --infile|-i   <input stream/file>\n",
        "    --padding|-p  [0+ digits of 0-padding]\n",
        "    --start|-s    [starting digit value; default \"1\"]\n",
        "    --help|-h\n",
        ;
}

if ( (! looks_like_number($pad) ) or ( $pad < 0 ) or ( $pad != int($pad) ) ) { 
    die "Digits of 0-padding must be non-negative integer, not: $pad\n";
}

if ( (! looks_like_number($i) ) or ( $i < 1 ) or ( $i != int($i) ) ) {
    die "Starting number must be positive integer, not: $i\n";
}

# Accept either a stream from '-' or a standard file.
my $INPUT_FILE;
if ($infile eq '-') {
    # Special case: get the stdin handle
    $INPUT_FILE = *STDIN{IO};
}
else {
    # Standard case: open the file
    open $INPUT_FILE, '<', $infile;
}

while (my $input = <$INPUT_FILE>) { 
    chomp $input;
    if ( $input !~ /\A > /xms ) { 
        print "$input\n";
    }
    elsif ( $input =~ /\A > (\S .*) \z/xms ) {
        $header = $1;

        if ( $pad >= 1 ) {
            my $format = '%0' . $pad . 'u';
            $index = sprintf "$format", $i;    # This lets me feed user-specified values into sprintf.
        }
        else {
            $index = $i;
        }

        print ">$index  $header\n";
        $i++;
    }
    else { 
        die "Can't parse input line: $input\n";
    }
}

close $INPUT_FILE;

