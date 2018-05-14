#!/usr/bin/env perl

# make_LAfix_script.pl -- auxillary Perl script for making serial LAfix commands in MARVEL assembly pipeline.
# Erich Schwarz <ems394@cornell.edu>, 3/5/2018.

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

# Use the first argment to get a database name:
my $database = q{};
$database    = shift @ARGV if @ARGV;

# Use the second argument to get a read length:
my $length = q{};
$length    = shift @ARGV if @ARGV;

# Treat all of the other arguments as input files:
my @infiles = ();
@infiles = @ARGV if @ARGV;

# Enforce the database name being sane (wordlike):
if ( $database !~ /\A \w+ \z/xms ) {
    warn "Format: make_LAfix_script.pl  [database name]  [read length, positive integer]  [1+ input *.las files]\n";
    die "Did not supply a clean, word-like database name (\"$database\")";
}

# Enforce the read length being a positive integer:
if ( (! looks_like_number($length) ) or ( $length < 1 ) or ( $length != int($length) ) ) {
    warn "Format: make_LAfix_script.pl  [database name]  [read length, positive integer]  [1+ input *.las files]\n";
    die "Did not supply a positive-integer read length (\"$length\")";
}

# Enforce there being at least one infile:
if (! @infiles) {
    warn "Format: make_LAfix_script.pl  [database name]  [read length, positive integer]  [1+ input *.las files]\n";
    die "No *.las input files!";
}

my %file2block = ();

# verify name format and map infiles to index numbers.
# sample infile name:  HAWAII.1.las
foreach my $infile (@infiles) {
    if ( $infile !~ /\A \w+ \. \d+ \. las \z/xms ) { 
        die "Infile name not parsable by make_LAfix_script.pl: $infile\n";
    }
    elsif ( $infile =~ /\A \w+ \. (\d+) \. las \z/xms ) {
        my $block = $1;
        $file2block{$infile} = $block;
    }
}

# sort infiles by true numbering, not ASCII numbering:
@infiles = sort { $file2block{$a} <=> $file2block{$b} } @infiles;

# finally, print out line-commands for each infile:
foreach my $infile (@infiles) {
    if ( $infile =~ /\A (\w+ \. (\d+)) \. las \z/xms ) {
        my $stem    = $1;
        my $block   = $2;
        my $outfile = "$stem.fixed.fasta";
        print "LAfix -x $length $database $infile $outfile ;\n"
    }
}

