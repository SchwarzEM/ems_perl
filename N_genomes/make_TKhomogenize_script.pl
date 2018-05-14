#!/usr/bin/env perl

# make_TKhomogenize_script.pl -- auxillary Perl script for making serial TKhomogenize commands in MARVEL assembly pipeline.
# Erich Schwarz <ems394@cornell.edu>, 3/8/2018.

use strict;
use warnings;
use autodie;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my @infiles  = ();
my $database = q{};
my $start    = 'repeats';
my $end      = 'hrepeats';
my $track    = 'trim';

my $help;

# TKhomogenize arguments that we are populating:
# usage: TKhomogenize [-m] [-ber n] [-iIt track] database input.las
# [we need to tell it the DATABASE]
# also:
# -i track  input interval track (repeats)
# -I track  output interval track (hrepeats)
# -t track  trim annotation track to be used (trim)

GetOptions ( 'infiles=s{,}' => \@infiles,
             'database=s'   => \$database,
             'start=s',     => \$start,
             'end=s'        => \$end,
             'track=s'      => \$track,
             'help'         => \$help,   );

if ( $help 
     or (! @infiles) 
     or ( $database !~ /\A \w+ \z/xms ) 
     or ( $start    !~ /\A \w+ \z/xms )
     or ( $end      !~ /\A \w+ \z/xms )
     or ( $track    !~ /\A \w+ \z/xms )
   ) { 
    die "Format: make_LArepeat_script.pl\n",
        "    --infiles|-i    <input files>\n",
        "    --database|-d   [mandatory DATABASE argument]\n",
        "    --start|-s      [supplies \"-i\" value; input interval track; default value is \"repeats\"]\n",
        "    --end|-e        [supplies \"-I\" value; output interval track; default value is \"hrepeats\"]\n",
        "    --track|-t      [supplies \"-t\" value; trim annotation track to be used; default value is \"trim]\n",
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
        # sample line:  TKhomogenize -i repeats -I hrepeats -t trim0 -b <block> HAWAII_FIX HAWAII_FIX.<block>.stitch.las ;
        print "TKhomogenize -i $start -I $end -t $track -b $block $database $infile ;\n";
    }
}

