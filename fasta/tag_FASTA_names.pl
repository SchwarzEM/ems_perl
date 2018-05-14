#!/usr/bin/env perl

# tag_FASTA_names.pl -- Erich Schwarz, <ems394@cornell.edu>, 7/20/2013.
# Purpose: append prefixes or suffixes (or both) onto header names of one or more FASTA files; can act as filter on piped stream.

use strict;
use warnings;
use Getopt::Long;

my @infiles = ();
my $prefix  = q{};
my $suffix  = q{};
my $help;

GetOptions ( 'infile=s{,}' => \@infiles,
             'prefix=s'    => \$prefix,
             'suffix=s'    => \$suffix,
             'help'        => \$help,   );

if ( ($help) or (! @infiles ) or ( (! $prefix ) and (! $suffix) ) ) { 
    die "Format: tag_FASTA_names.pl --infile|-i <input files/stream> --prefix|-p [prefix] --suffix|-s [suffix] --help|-h [type this message]\n";
}

if ( $prefix and ( $prefix !~ / \A \w[\w\.\-\|_]* \z /xms ) ) { 
    die "Cannot use prefix $prefix!\n";
}
if ( $suffix and ( $suffix !~ / \A \w[\w\.\-\|_]* \z /xms ) ) { 
    die "Cannot use suffix $suffix!\n";
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
        if ($input !~ / \A > \S /xms) {
            print "$input\n";
        }
        elsif ($input =~ / \A > (\S+) (.*?) \z /xms) {
            my $input1 = $1;
            my $input2 = $2;
            print ">";
            print $prefix if $prefix;
            print $input1;
            print $suffix if $suffix;
            print $input2;
            print "\n";
        }
        else { 
            die "Can't parse input line: $input!\n";
        }
    }
    close $INPUT_FILE or die "Can't close filehandle to input file $infile. $!\n";
}
