#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $subdir_list = q{};
$subdir_list    = $ARGV[0] if $ARGV[0];

my $dir = '/ocean/projects/mcb190015p/schwarze/nippo_steiner/nippo/genbank/uploads/ucr/';

if (! $subdir_list) {
    die "Format: make_job_R1_listing_2023.09.14.01.pl [dir list] > [script for listing R1 read sets by dir]\n";
}

open my $SUBDIRS, '<', $subdir_list;
while ( my $subdir = <$SUBDIRS> ) {
    chomp $subdir;
    my $outfile = "$subdir.list";
    $outfile    = safename($outfile);
    print "ls /ocean/projects/mcb190015p/schwarze/nippo_steiner/nippo/genbank/uploads/ucr/$subdir/*.R1.fastq.gz > $outfile ;\n";
}
close $SUBDIRS;


sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}

