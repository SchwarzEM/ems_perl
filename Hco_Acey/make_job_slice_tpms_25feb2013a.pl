#!/usr/bin/env perl

# make_job_slice_tpms_25feb2013a.pl -- Erich Schwarz <ems394@cornell.edu>, 2/25/2013.
# Purpose: given a list of RSEM 1.2.0 gene results files, make a script that will efficiently slice them; don't put hardwired Cel_ prefix on them.

use strict;
use warnings;
use File::Basename;

while (my $input = <>) {
    chomp $input;
    my $stem = basename $input;
    $stem =~ s/\.genes\.results//;
    print "    /mnt/home/emsch/perl.svn/trunk/afd2solexa_data/extract_rsem_slices.pl -p -o $stem -i $input -s 9 ;\n";
}

