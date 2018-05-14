#!/usr/bin/env perl

# md5sum_recheck.pl -- Erich Schwarz <ems394@cornell.edu>, 3/26/2013.
# Purpose: given many MD5sum check files, verify that files with their names have the same md5sum score.  Pipe contents in.

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \s+ (\S+) \s* \z/xms ) { 
        my $putative_md5sum = $1;
        my $file            = $2;
        if (-r $file) { 
            my $result = `md5sum $file | grep $putative_md5sum`;
            chomp $result;
            if ( $result !~ /\S/xms ) { 
                print "Failed to verify md5sum of $file!\n";
            }
            else { 
                print "Verified md5sum of $file: $result\n";
            }
        }
    }
}

