#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while ( my $bedgraph = <> ) {
    chomp $bedgraph;
    if ( $bedgraph =~ /\A (\S+) \. bedgraph \z/xms ) {
        my $stem  = $1;
        my $sizes = q{};

        my $bw    = "$stem.bw";
        $bw       = safename($bw);

        if ( $stem =~ /\A Necator /xms ) {
            $sizes = '/ocean/projects/mcb190015p/shared/schwarze/Necator/2024.05.11/aroian/Necator_2022.05.29.02.chrs_only.chrom.sizes';
        }
        elsif ( $stem =~ /\A Ilik2 /xms ) {
            $sizes = '/ocean/projects/mcb190015p/shared/schwarze/Necator/2024.05.11/ilik2/Ilik2_2024.02.23.chrs_only.dna.chrom.sizes';
        }
        else {
            die "Can't assign stem $stem to *.chrom.sizes files\n";
        }

        print '$PROJECT/src/kent_v362_linux.x86_64/bedGraphToBigWig ' . "$bedgraph $sizes $bw ;\n";
    }
    else {
        die "Can't parse input BedGraph file name: $bedgraph\n";
    }
}

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

