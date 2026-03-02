#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;

my $list    = q{};
my $suffix  = q{};

$list       = $ARGV[0] if $ARGV[0];
$suffix     = $ARGV[1] if $ARGV[1];

if ( (! $list ) or (! $suffix ) ) {
    die "Format: make_mosdepths_02mar2026.pl [input BAM list] [suffix] > [mosdepth commands]\n";
}

open my $LIST, '<', $list;
while ( my $infile = <$LIST> ) {
    chomp $infile;

    if (! -e $infile) {
        die "From infile list $list, cannot find infile: $infile\n";
    }

    my $basename = basename($infile);
    my $stem = q{};

    if ( $basename =~ /\A (\S+) \. bam \z/xms ) {
        $stem = $1;
        my $output = "$stem.$suffix";
        print "mosdepth --threads 16 --no-per-base --by 10000 $output $infile ;\n";
    }
    else {
        die "Cannot parse stem from infile basename $basename\n";
    }
}
close $LIST;
