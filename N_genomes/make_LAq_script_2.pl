#!/usr/bin/env perl

# make_LAq_script_2.pl -- auxillary Perl script for making serial LAq commands in MARVEL assembly pipeline.
# Erich Schwarz <ems394@cornell.edu>, 3/7/2018.

use strict;
use warnings;
use autodie;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my $database = q{};
my $track    = q{};
my @infiles  = ();

my $help;

GetOptions ( 'infiles=s{,}' => \@infiles,
             'database=s'   => \$database,
             'track=s'      => \$track,
             'help'         => \$help,   );

if ( $help 
     or (! @infiles) 
     or ( $database !~ /\A \w+ \z/xms ) 
     or ( $track and ( $track !~ /\A \w+ \z/xms ) ) 
   ) { 
    die "Format: count_fasta_residues.pl\n",
        "    --infiles|-i    <input files>\n",
        "    --database|-d   [mandatory DATABASE argument]\n",
        "    --track|-t      [optional track argument]\n",
        "    --help|-h       [print this message]\n",
        ;
}

my %file2block = ();

# verify name format and map infiles to index numbers.
# sample infile name:  HAWAII_FIX.<block>.stitch.las, e.g., HAWAII_FIX.1.stitch.las
# but also design this so that it should work with HAWAII_FIX.<block>.las!

foreach my $infile (@infiles) {
    if ( $infile !~ /\A \w+ \. \d+ (?: \. stitch){0,1} \. las \z/xms ) { 
        die "Infile name not parsable by make_LAq_script_2.pl: $infile\n";
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
        # sample line:  LAq -b <block> HAWAII_FIX HAWAII_FIX.<block>.stitch.las ;
        print "LAq -b $block";
        if ($track) {
            print " -T $track";
        }
        print " $database $infile ;\n";
    }
}

