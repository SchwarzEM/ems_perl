#!/usr/bin/env perl

# make_LAq_script_3.pl -- auxillary Perl script for making serial LAq commands in MARVEL assembly pipeline.
# Erich Schwarz <ems394@cornell.edu>, 3/8/2018.

use strict;
use warnings;
use autodie;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my @infiles  = ();
my $database = q{};
my $start    = 'trim';
my $end      = 'trim';

my $help;

# LAq arguments that we are working with:
# usage: [-u] [-b n] [-d n] [-s n] [-S n] [-t track] [-T track] [-q track] [-Q track] database input.las
# -u        update existing trim track
# -b n      block number
# -t track  input trim track in -u mode (default trim)
# -T track  output trim track (default trim)

GetOptions ( 'infiles=s{,}' => \@infiles,
             'database=s'   => \$database,
             'start=s'      => \$start,
             'end=s'        => \$end,
             'help'         => \$help,   );

if ( $help 
     or (! @infiles) 
     or ( $database !~ /\A \w+ \z/xms )
     or ( $start    !~ /\A \w+ \z/xms )
     or ( $end      !~ /\A \w+ \z/xms )
   ) { 
    die "Format: count_fasta_residues.pl\n",
        "    --infiles|-i    <input files>\n",
        "    --database|-d   [mandatory DATABASE argument]\n",
        "    --start|-s      [supplies \"-t\" value; input trim track in -u (update) mode; default \"trim\"]\n",
        "    --end|-e        [supplies \"-T\" value; output trim track; default \"trim\"]\n",
        "    --help|-h       [print this message]\n",
        ;
}

my %file2block = ();

# verify name format and map infiles to index numbers.
# sample infile name:  HAWAII_FIX.<block>.gap.las, e.g., HAWAII_FIX.1.gap.las

foreach my $infile (@infiles) {
    if ( $infile !~ /\A \w+ \. \d+ \. gap \. las \z/xms ) { 
        die "Infile name not parsable by make_make_LAq_script_3.pl : $infile\n";
    }
    elsif ( $infile =~ /\A \w+ \. (\d+) \. gap \. las \z/xms ) {
        my $block = $1;
        $file2block{$infile} = $block;
    }
}

# sort infiles by true numbering, not ASCII numbering:
@infiles = sort { $file2block{$a} <=> $file2block{$b} } @infiles;

# finally, print out line-commands for each infile:
foreach my $infile (@infiles) {
    if ( $infile =~ /\A \w+ \. (\d+) \. gap \. las \z/xms ) {
        my $block   = $1;
        # sample line:  LAq -u -t trim0 -T trim1 -b <block> HAWAII_FIX HAWAII_FIX.<block>.gap.las ;
        print "LAq -u -t $start -T $end -b $block $database $infile ;\n";
    }
}

