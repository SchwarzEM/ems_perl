#!/usr/bin/perl

# scriptify_integr8_picks.pl, Erich Schwarz <emsch@its.caltech.edu>, 2/24/06.
# Purpose: convert selected taxon descriptions to script commands, e.g. --
# 
# 100     204722  Brucella suis   B_suis  325
#
# to
#
# wget ftp://ftp.ebi.ac.uk/pub/databases/integr8/uniprot/proteomes/100.B_suis.dat.gz ;
# gunzip 100.B_suis.dat.gz ;

use strict;
use warnings;

while (<>) { 
    chomp(my $input = $_);
    unless ($input =~ /^# /) { 
        $input =~ m/(\d+)\t\d+\t[^\t]+\t(\S+)/;
        my $arb_no   = $1;
        my $spec_tag = $2;
        print 'wget ftp://ftp.ebi.ac.uk/pub/databases/integr8/uniprot/proteomes/';
        print "$arb_no" . '.' . "$spec_tag" . '.dat.gz ;' . "\n";
        print "gunzip $arb_no" . '.' . "$spec_tag" . '.dat.gz ;' . "\n";
    }
}

