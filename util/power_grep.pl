#!/usr/bin/env perl

# power_grep.pl -- Erich Schwarz <ems394@cornell.edu>, 2/21/2013.
# Purpose: given a file of patterns and an incoming data stream, filter the stream +/- the patterns.

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw(uniq);

my @input_files  = ();
my $pattern_file = q{};
my @patterns     = ();
my $v;
my $help;

GetOptions ( 'input_files=s{,}'  => \@input_files,
             'v'                 => \$v,
             'pattern_file:s'    => \$pattern_file,
             'help'              => \$help,   );

if ( ($help) or (! @input_files) or (! $pattern_file ) ) { 
    die "Format: power_grep.pl",
        " --input_files|-i [input stream/files] --pattern_file|-p [file of pattern (regex) lines] --v|-v [pass NONmatches]\n",
        ;
}

# Upload pattern lines, and remove redundant ones.
open my $PATT, '<', $pattern_file or die "Can't open pattern list file $pattern_file: $!";
while (my $input = <$PATT>) { 
    chomp $input;
    $input = qr/$input/;
    push @patterns, $input;
}
close $PATT or die "Can't close filehandle to pattern list file $pattern_file: $!";
@patterns = uniq @patterns;

my $INPUT_FILE;

foreach my $infile (@input_files) { 
    # Allow multifile and/or streamed input.
    if ($infile eq '-') {
        # Special case: get the stdin handle
        $INPUT_FILE = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $INPUT_FILE, '<', $infile or die "Can't open input file $infile. $!\n";
    }

LOOP: while (my $input = <$INPUT_FILE>) { 
        my $to_print = 0;
        foreach my $pattern (@patterns) { 
            if ( ( $input =~ /$pattern/ ) and (! $v) ) { 
                print $input;
                next LOOP;
            }
        }
        print $input if ($v);
    }
    close $INPUT_FILE or die "Can't close filehandle to input file $infile: $!\n";
}

