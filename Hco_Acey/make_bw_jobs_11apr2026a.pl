#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $bedgraphs  = q{};
my $sizes_file = q{};

$bedgraphs  = $ARGV[0] if $ARGV[0];
$sizes_file = $ARGV[1] if $ARGV[1];

if ( (! $bedgraphs ) or (! $sizes_file ) ) {
    die "Format: make_bw_jobs_11apr2026a.pl [bedgraph file list] [sizes file] > [bedGraphToBigWig commands]\n";
}

open my $BEDGRAPHS, '<', $bedgraphs;
while ( my $bedgraph = <$BEDGRAPHS> ) {
    chomp $bedgraph;
    if ( $bedgraph =~ /\A (\S+) \. bedgraph \z/xms ) {
        my $stem  = $1;
        my $bw    = "$stem.bw";
        $bw       = safename($bw);
        print '$PROJECT/src/kent_v362_linux.x86_64/bedGraphToBigWig ' . "$bedgraph $sizes_file $bw ;\n";
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

