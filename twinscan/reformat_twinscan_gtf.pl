#!/usr/bin/perl

# reformat_twinscan_gtf.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/8/2007.
# Purpose: reformat Twinscan *.gtf files to be read as GFF by Bio::DB::GFF.  See comments below.

use strict;
use warnings;

my ($input, $seq_length, $seq_name, $watching);
$watching = 0;

$^I = ".bak";

while ($input = <>) { 
    chomp $input;
    if ($input =~ /^#/) { 
        $watching = 1;
    }
    if ($input =~ /^# Target Sequence: >(\S+)/) { 
        $seq_name = $1;
    }
    if ($input =~ /^# Target Sequence Read... (\d+)bp/) { 
        $seq_length = $1;
    }
    if (($input !~ /^#/) and ($watching)) { 
        print "$seq_name\tassembly\tcontig\t";
        print "1\t$seq_length\t.\t+\t.\t";
        print "Sequence $seq_name\n";
        $watching = 0;
    }
    print "$input\n";
}


# 
# The reason for all this is in the Bio::DB::GFF documentation:
# 
# """
#
# ... each annotation in a GFF file refers to a reference sequence.  It is
# important that each reference sequence also be identified by a line in the
# GFF file.  This allows the Bio::DB::GFF module to determine the length and
# class of the reference sequence, and makes it possible to do relative
# arithmetic.
# 
# For example, if "Chr1" is used as a reference sequence, then it should
# have an entry in the GFF file similar to this one:
# 
# Chr1 assembly chromosome 1 14972282 . + . Sequence Chr1
# 
# This indicates that the reference sequence named "Chr1" has length
# 14972282 bp, method "chromosome" and source "assembly".
# 
# In addition, as indicated by the group field, Chr1 has class "Sequence"
# and name "Chr1".
# 
# """
# 

