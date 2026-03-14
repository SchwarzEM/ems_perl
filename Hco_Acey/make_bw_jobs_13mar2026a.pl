#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while ( my $bed = <> ) {
    chomp $bed;
    if ( $bed =~ /\A (\S+) \. bed \z/xms ) {
        my $stem  = $1;
        my $bgrh  = "$stem.bedgraph";
        my $s_bg  = "$stem.sorted.bedgraph";
        my $bw    = "$stem.bw";
        my $fai   = q{};
        my $sizes = q{};

        if ( $stem =~ /\A Necator /xms ) {
            $fai   = '/ocean/projects/mcb190015p/shared/schwarze/Necator/2024.05.11/aroian/Necator_2022.05.29.02.chrs_only.dna.fa.fai';
            $sizes = '/ocean/projects/mcb190015p/shared/schwarze/Necator/2024.05.11/aroian/Necator_2022.05.29.02.chrs_only.chrom.sizes';
        }
        elsif ( $stem =~ /\A Ilik2 /xms ) {
            $fai   = '/ocean/projects/mcb190015p/shared/schwarze/Necator/2024.05.11/ilik2/Ilik2_2024.02.23.chrs_only.dna.fa.fai';
            $sizes = '/ocean/projects/mcb190015p/shared/schwarze/Necator/2024.05.11/ilik2/Ilik2_2024.02.23.chrs_only.dna.chrom.sizes';
        }
        else {
            die "Can't assign stem $stem to *.fai and *.chrom.sizes files\n";
        }

        print "mamba activate bedtools_2.30.0 ;\n";
        print "bedtools genomecov -bga -i $bed -g $fai > $bgrh ;\n";
        print "mamba deactivate ;\n";
        print "sort -k1,1 -k2,2n $bgrh > $s_bg ;\n";
        print '$PROJECT/src/kent_v362_linux.x86_64/bedGraphToBigWig ' . "$s_bg $sizes $bw ;\n";
    }
    else {
        die "Can't parse input BED file name: $bed\n";
    }
}

