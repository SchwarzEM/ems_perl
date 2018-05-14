#!/usr/bin/env perl

# make_LArepeat_script.pl -- auxillary Perl script for making serial LArepeat commands in MARVEL assembly pipeline.
# Erich Schwarz <ems394@cornell.edu>, 3/8/2018.

use strict;
use warnings;
use autodie;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my @infiles  = ();
my $database = q{};
my $coverage = q{};
my $lower    = 1.7;
my $upper    = 2.0;

my $help;

# LArepeat arguments that we are populating:
# usage: [-h f] [-l f] [-t track] [-b n] [-c n] [-m n] [-n n] [-o n] DATABASE input.las
# [we need to tell it the DATABASE]
# also:
# -c n  expected coverage of the dataset. -1 auto-detect. (-1)
# -l f  below which multiple of the expected coverage a repeat ends (1.7)
# -h f  above which multiple of the expected coverage the start of a repeat is reported (2.0)

GetOptions ( 'infiles=s{,}' => \@infiles,
             'database=s'   => \$database,
             'coverage=f',  => \$coverage,
             'lower=f'      => \$lower,
             'upper=f'      => \$upper,
             'help'         => \$help,   );

if ( $help 
     or (! @infiles) 
     or ( $database !~ /\A \w+ \z/xms ) 
     or (! looks_like_number($coverage) ) 
     or (! looks_like_number($lower) ) 
     or (! looks_like_number($upper) )
   ) { 
    die "Format: make_LArepeat_script.pl\n",
        "    --infiles|-i    <input files>\n",
        "    --database|-d   [mandatory DATABASE argument]\n",
        "    --coverage|-c   [coverage of genome (note that for size-filtered reads, *not* defined by bulk read nt count/genome size)]\n",
        "    --lower|-l      [below which multiple of the expected coverage a repeat ends; program default is 1.7]\n",
        "    --upper|-u      [supplies \"-h\" value; above which multiple of the expected coverage the start of a repeat is reported; default is 2.0]\n",
        "    --help|-h       [print this message]\n",
        ;
}

my %file2block = ();

# verify name format and map infiles to index numbers.
# sample infile name:  HAWAII_FIX.<block>.stitch.las, e.g., HAWAII_FIX.1.stitch.las
# but also design this so that it should work with HAWAII_FIX.<block>.las!

foreach my $infile (@infiles) {
    if ( $infile !~ /\A \w+ \. \d+ (?: \. stitch){0,1} \. las \z/xms ) { 
        die "Infile name not parsable by make_LArepeat_script.pl: $infile\n";
    }
    elsif ( $infile =~ /\A \w+ \. (\d+) (?: \. stitch){0,1} \. las \z/xms ) {
        my $block = $1;
        $file2block{$infile} = $block;
    }
}

# sort infiles by true numbering, not ASCII numbering:
@infiles = sort { $file2block{$a} <=> $file2block{$b} } @infiles;

# finally, print out line-commands for each infile:
foreach my $infile (@infiles) {
    if ( $infile =~ /\A \w+ \. (\d+) (?: \. stitch){0,1} \. las \z/xms ) {
        my $block   = $1;
        # sample line:  LArepeat -c <coverage> -l 1.5 -h 2.0 -b <block> HAWAII_FIX HAWAII_FIX.<block>.stitch.las ;
        print "LArepeat -c $coverage -l $lower -h $upper -b $block $database $infile ;\n";
    }
}

