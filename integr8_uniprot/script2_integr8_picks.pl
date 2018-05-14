#!/usr/bin/perl

# script1_integr8_picks.pl, Erich Schwarz <emsch@its.caltech.edu>, 2/24/06.
# Purpose: convert selected taxon descriptions to script commands, e.g. --
# 
# 100     204722  Brucella suis   B_suis  325
#
# to
#
# mv -i   DRD_integr8_files/proteomes/100.B_suis.dat.gz   DRD_integr8_files/preferred_proteomes ;

use strict;
use warnings;

print 'rm    -rf DRD_integr8_files/preferred_proteomes ; ';
print "\n";
print 'mkdir -p  DRD_integr8_files/preferred_proteomes ; ';
print "\n\n";

while (<>) { 
    chomp(my $input = $_);
    unless ($input =~ /^# /) { 
        $input =~ m/(\d+)\t\d+\t[^\t]+\t(\S+)/;
        my $arb_no   = $1;
        my $spec_tag = $2;
        print 'mv -i   DRD_integr8_files/proteomes/';
        print "$arb_no";
        print '.';
        print "$spec_tag";
        print '.dat.gz   DRD_integr8_files/preferred_proteomes';
        print " ; \n";
    }
}

