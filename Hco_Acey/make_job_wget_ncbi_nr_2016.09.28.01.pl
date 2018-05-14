#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = '#!/bin/bash' . "\n";
print "$header\n";

foreach my $i (0..55) {
    my $j = sprintf "%02i", $i;

    print "wget -q ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.$j.tar.gz ;\n";
    print "wget -q ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.$j.tar.gz.md5 ;\n";
    print "\n";
}

