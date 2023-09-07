#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A \S+ \z/xms ) {
        print "mkdir -p /ocean/projects/mcb190015p/shared/schwarze/Nippostronglyus_project/fastq/$input ;\n";
        print "rsync -av ",
              "/ocean/projects/mcb190015p/isaryhia/Nippostronglyus_project/fastq/$input/* ",
              "/ocean/projects/mcb190015p/shared/schwarze/Nippostronglyus_project/fastq/$input ;\n",
              ;
    }
    else {
        die "Cannot parse input: $input\n";
    }
}

