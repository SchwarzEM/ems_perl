#!/usr/bin/env perl

# filter_func_refinement.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/3/2010.
# Purpose: select only parts of refined FUNC output which have an actual result;

use strict;
use warnings;
use Getopt::Long;

my @input_files = ();

my $upreg;
my $downreg;
my $help;

GetOptions ( 'input=s{,}' => \@input_files,
             'upreg'      => \$upreg,
             'downreg'    => \$downreg,
             'help'       => \$help, );

if ( $help or (! @input_files ) or ( $upreg and $downreg ) ) { 
    die "Format: filter_func_refinement.pl\n",
        "    --input      [1+ files, or '-' for input stream]\n",
        "    --upreg|-u   [upregulated only]\n",
        "    --downreg|-d [downregulated genes only]\n",
        "    --help|-h    [print this message]\n",
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
        if ( $input =~ / \A 
                         [^\t]+
                         \t [^\t]+
                         \t GO:\d+ 
                         \t \+ 
                         \t [^\t]* 
                         \t [^\t]* 
                         \t ([^\t]+) 
                         \t ([^\t]+) 
                         \t 
                         \z /xms ) {
            my $downreg_prob = $1;
            my $upreg_prob   = $2;
            if ( (! $upreg and ( $upreg_prob > 0.5 ) ) or (! $downreg and ( $downreg_prob > 0.5 ) ) ) { 
                print "$input\n";
            }
        }
        elsif ( $input =~ /\A 
                           root_node_name
                           \t node_name
                           \t node_id 
                           \t sign\?
                           \t raw_p_low_ranks 
                           \t raw_p_high_ranks 
                           \t p_low_ranks_after_refinement 
                           \t p_high_ranks_after_refinement 
                           \z /xms ) { 
            print "$input\n";
        }
        elsif ( $input !~ / \A 
                            [^\t]+
                            \t [^\t]+
                            \t GO:\d+
                            \t \-
                            \t [^\t]*
                            \t [^\t]*
                            \t [^\t]+
                            \t [^\t]+
                            \t
                            \z /xms ) {
            die "Can't parse input line: $input!\n";
        }
    }
    close $INPUT_FILE or die "Can't close filehandle to input file $infile. $!\n";
}

