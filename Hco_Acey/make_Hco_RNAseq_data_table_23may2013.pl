#!/usr/bin/env perl

use strict;
use warnings;

my %long = ( 'Egg' => 'egg stage',
             'L1' => 'first-stage larvae (L1)',
             'L2-1' => 'second-stage larvae (L2), biological replicate 1',
             'L2-2' => 'second-stage larvae (L2), biological replicate 2',
             'L4.female' => 'female fourth-stage larvae (L4)',
             'L4.male' => 'male fourth-stage larvae (L4)',
             'Adult.female' => 'female adult stage',
             'Adult.male' => 'male adult stage', );

while (my $input = <>) { 
    chomp $input;
    if ( $input !~ /\A \S+ \s+ Hco_\S+\.lib_\d{5}\.\S+\.\d\.fq\.gz \z/xms ) { 
        die "Can't parse input: $input\n";
    }
    if ( $input =~ /\A (\S+) \s+ (Hco_(\S+)\.lib_(\d{5})\.(\S+)\.(\d)\.fq\.gz) \z/xms ) {
        my $md5  = $1;
        my $file = $2;
        my $type = $3;
        my $lib  = $4;
        my $cell = $5;
        my $lane = $6;

        if (! exists $long{$type}) { 
            die "Can't expand type: $type\n";
        }

        my $printed_type = $type;
        $printed_type    =~ s/\./ /g;

        print "File: $file\n";
        print "md5sum: $md5\n";
        print '#', "\n";
        print "Alias: Haemonchus contortus (Haecon-5), $printed_type\n";
        print "Title: RNA-seq from Haemonchus contortus (Haecon-5), $long{$type}\n";
        print "NCBI Taxon ID: 6289\n";
        print "Description: RNA-seq reads from Haemonchus contortus (Haecon5), $long{$type};\n";
        print "Jacobs Genome Center Illumina sequencing library $lib, flowcell $cell, lane $lane.\n";
        print "\n";
    }
}

