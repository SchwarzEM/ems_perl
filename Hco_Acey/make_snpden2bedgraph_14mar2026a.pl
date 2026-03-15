#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

my $list = q{};
my $bin  = q{};
my $chr  = q{};

$list = $ARGV[0] if $ARGV[0];
$bin  = $ARGV[1] if $ARGV[1];
$chr  = $ARGV[2] if $ARGV[2];

if ( (! $list ) or (! $bin ) or (! looks_like_number($bin) ) or (! $chr ) ) {
    die "Format: make_snpden2bed_13mar2026a.pl [list of input snpden files] [bin size] [chr size file] > [snpden2bedgraph_14mar2026.pl commands]\n";
}

open my $LIST, '<', $list;
while ( my $infile = <$LIST> ) {
    chomp $infile;
    if ( $infile =~ /\A (\S+) \. snpden \z/xms ) {
        my $stem = $1;
        my $rawfile = "$stem.orig.bedgraph";
        my $outfile = "$stem.bedgraph";
        my $err     = "$stem.err";
        print '$PROJECT/ems_perl/Hco_Acey/snpden2bedgraph_14mar2026.pl ';
        print "$infile $bin $chr 1>$rawfile 2>$err ;\n";
        print "sort -k1,1 -k2,2n $rawfile > $outfile ;\n";
    }
    else {
        die "Cannot parse infile: $infile\n";
    }
}
close $LIST;

