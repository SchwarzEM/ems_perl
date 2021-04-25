#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    # sample input: /ocean/projects/mcb190015p/schwarze/Acey/scripts/job_acey_salmon_2019.05.28.01.sh
    if ( $input =~ /\A\S+\/([^\/\s]+)\z/xms ) { 
        my $output = $1;
        $output =~ s/2019\.05\.28/2021.04.24/g;
        print "   cat $input | ./update_salmon_24apr2021.pl > $output ;\n";
    }
    else {
        die "Cannot parse: $input\n"
    }
}
