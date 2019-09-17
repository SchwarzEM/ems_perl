#!/usr/bin/env perl

# rescue_horrible_fasta_names.pl -- Erich Schwarz <ems394@cornell.edu>, 9/17/2019.
# Purpose: take horribly misheadered FASTAs (e.g., from 'polishing' by a certain famous genome assembler) and make them sane again.

use strict;
use warnings;
use autodie;

use Getopt::Long;

my $i            = 0;
my $prefix       = q{};
my $infile       = q{};
my $seqname      = q{};
my $help;

GetOptions ( 'infile=s'   => \$infile,
             'prefix:s'   => \$prefix,
             'help'       => \$help,     );

# Have human-readable default value:
$prefix   ||= 'Scaffold.';

if ( $help or (! $infile ) ) { 
    die "Format: rescue_horrible_fasta_names.pl\n",
        "    --infile|-i    <input stream/files>\n",
        "    --prefix|-p    [optional prefix; default is \"Scaffold_\"]\n",
        "    --help|-h\n",
}

# Accept either a stream from '-' or a standard file.
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
    if ( $input =~ /\A > (.*) \z/xms ) { 
        my $header = $1;
        $i++;
        $seqname = $prefix . $i;
        print ">$seqname  $header\n";
    }
    else {
        print "$input\n";
    }
}
close $INPUT_FILE;

