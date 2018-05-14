#!/usr/bin/env perl

# make_LAgap_script.pl -- auxillary Perl script for making serial LAgap commands in MARVEL assembly pipeline.
# Erich Schwarz <ems394@cornell.edu>, 3/8/2018.

use strict;
use warnings;
use autodie;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my @infiles  = ();
my $database = q{};
my $stitch   = 0;
my $track    = q{};

my $help;

# LAgap arguments that we are populating:
# usage: LAgap [-pL] [-s n] [-e track] [-t track] database input.las output.las
# [we need to tell it the DATABASE]
# also:
# -s n  don't count alignments that would be stitchable with a maximum distance of n (0)
# -t track  trim overlaps before gap detection

GetOptions ( 'infiles=s{,}' => \@infiles,
             'database=s'   => \$database,
             'stitch=i',    => \$stitch,
             'track=s'      => \$track,
             'help'         => \$help,   );

if ( $help 
     or (! @infiles) 
     or ( $database !~ /\A \w+ \z/xms ) 
     or (! looks_like_number($stitch) )  
     or ( $stitch < 0                 )
     or ( $stitch != int($stitch)     ) 
     or ( $track    !~ /\A \w+ \z/xms )
   ) { 
    die "Format: make_LArepeat_script.pl\n",
        "    --infiles|-i    <input files>\n",
        "    --database|-d   [mandatory DATABASE argument]\n",
        "    --stitch|-s     [don't count alignments that would be stitchable with a maximum distance of n; default value is 0]\n",
        "    --track|-t      [supplies track value; trim overlaps before gap detection]\n",
        "    --help|-h       [print this message]\n",
        ;
}

my %file2block = ();

# verify name format and map infiles to index numbers.
# sample infile name:  HAWAII_FIX.<block>.stitch.las, e.g., HAWAII_FIX.1.stitch.las

foreach my $infile (@infiles) {
    if ( $infile !~ /\A \w+ \. \d+ \. stitch \. las \z/xms ) { 
        die "Infile name not parsable by make_LArepeat_script.pl: $infile\n";
    }
    elsif ( $infile =~ /\A \w+ \. (\d+) \. stitch \. las \z/xms ) {
        my $block = $1;
        $file2block{$infile} = $block;
    }
}

# sort infiles by true numbering, not ASCII numbering:
@infiles = sort { $file2block{$a} <=> $file2block{$b} } @infiles;

# finally, print out line-commands for each infile:
foreach my $infile (@infiles) {
    if ( $infile =~ /\A (\w+ \. (\d+)) \. stitch \. las \z/xms ) {
        my $stem    = $1;
        my $block   = $2;
        my $outfile = "$stem.gap.las";
        # sample line:  LAgap -s 300 -t trim0 HAWAII_FIX HAWAII_FIX.<block>.stitch.las HAWAII_FIX.<block>.gap.las ;
        print "LAgap -s $stitch -t $track $database $infile $outfile ;\n";
    }
}

