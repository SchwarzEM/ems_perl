#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %seen = ();

while (my $input = <>) {
    # Sample input lines:
    # Strain	Strain Number	Condition	Treatment	Plate	Replicate	File	Folder
    # T727	TWF154	T727-10hr	10 hr after N2 added	LNM	1	32-T727-10-1_S5_R1_001.fastq.gz	N414-20170619-TruSeq_RNA

    chomp $input;
    if ( $input =~ /\A T727 \t .+ \t (\S+\.fastq\.gz) \t (\S+) \z/xms ) { 
        my $basename = $1;
        my $subdirectory = $2;

        if ( $seen{$basename} ) {
            die "Redundant input file: $basename\n";
        }

        my $path_to_fullname = '/pylon5/mc5fpnp/gvidald/N414_RNA_Data/' . $subdirectory . '/*/*/' . $basename;
        my $fullname = `ls $path_to_fullname`;
        chomp $fullname;
        print "$fullname\n";
    }
    elsif ( $input !~ /\A Strain /xms ) {
        die "Cannot parse: $input\n";
    }
}

