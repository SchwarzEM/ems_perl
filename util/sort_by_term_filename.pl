#!/usr/bin/env perl

# sort_by_term_filename.pl -- Erich Schwarz <ems394@cornell.edu>, 1/24/2016.
# Purpose: sort Unix etc. file lists by the file names alone, not their directories.
# Note: *almost* simple enough to be an inline Perl command, but complicated and fragile enough to be a small script instead.

use strict;
use warnings;
use autodie;

use File::Basename;

my @input_files = <>;

# Schwartzian transform; ref.: https://en.wikipedia.org/wiki/Schwartzian_transform
my @output_files = map  { $_->[0] }
                   sort { $a->[1] cmp $b->[1] }
                   map  { [$_, basename($_)] }
                   @input_files ;

print @output_files;

