#!/usr/bin/perl

# script1_integr8_picks.pl, Erich Schwarz <emsch@its.caltech.edu>, 2/24/06.
# Purpose: convert selected taxon descriptions to script commands, e.g. --
# 
# 100     204722  Brucella suis   B_suis  325
#
# to
#
# mv DRD_integr8_files/proteomes/100.B_suis.dat ../unused_proteomes ;

use strict;
use warnings;

print "mkdir -p DRD_integr8_files/unused_proteomes ; \n\n";

while (<>) { 
    chomp(my $input = $_);
    if ($input =~ m/^#\s+(\d+)\t\d+\t[^\t]+\t(\S+)/) { 
        my $arb_no   = $1;
        my $spec_tag = $2;
        print 'mv -i   DRD_integr8_files/proteomes/';
        print "$arb_no";
        print '.';
        print "$spec_tag";
        print '.dat.gz   DRD_integr8_files/unused_proteomes';
        print " ; \n";
    }
}

