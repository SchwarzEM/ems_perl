#!/usr/bin/env perl

# make_LAq_script.pl -- auxillary Perl script for making serial LAq commands in MARVEL assembly pipeline.
# Erich Schwarz <ems394@cornell.edu>, 3/4/2018.

use strict;
use warnings;
use autodie;

my @infiles = @ARGV;

my %file2block = ();

# verify name format and map infiles to index numbers:
foreach my $infile (@infiles) {
    if ( $infile !~ /\A \w+ \. \d+ \. las \z/xms ) { 
        die "Infile name not parsable by make_LAq_script.pl: $infile\n";
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
    if ( $infile =~ /\A \w+ \. (\d+) \. las \z/xms ) {
        my $block = $1;
        print "LAq -b $block HAWAII $infile ;\n"
    }
}
