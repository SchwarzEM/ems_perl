#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Cwd;

my $dir = getcwd;

my $header = '#!/bin/bash' . "\n\n" . "cd $dir;\n\n";

while (my $input = <>) { 
    chomp $input;
    my $basename = basename($input);
    print $header if $header;
    $header = q{};
    print "../tbl2asn -p . -t ../Acey_15jan2014.sbt",
          " -i $input",
          ' -j "[organism=Ancylostoma ceylanicum] [tech=wgs] [strain=HY135]"',
          " -l paired-ends -a r1k -c f -V vb -M n -X E",
          " -Z discrep_", $basename, ".txt ;",
          "\n",
          ;
}

print "\n";

# Censor this, to avoid confusing NCBI.
# print "    e_ping -p done_make_mass_tbl2asn_23feb2014 ;\n\n";


