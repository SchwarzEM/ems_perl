#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

my $list = q{};
my $bin  = q{};

$list = $ARGV[0] if $ARGV[0];
$bin  = $ARGV[1] if $ARGV[1];

if ( (! $list ) or (! $bin ) or (! looks_like_number($bin) ) ) {
    die "Format: make_snpden2bed_13mar2026a.pl [list of input snpden files] [bin size] > [snpden2bed_13mar2026.pl commands]\n";
}

open my $LIST, '<', $list;
while ( my $infile = <$LIST> ) {
    chomp $infile;
    if ( $infile =~ /\A (\S+) \. snpden \z/xms ) {
        my $stem = $1;
        my $outfile = "$stem.bed";
        my $err     = "$stem.err";
        print '$PROJECT/ems_perl/Hco_Acey/snpden2bed_13mar2026.pl ';
        print "$infile $bin 1>$outfile 2>$err ;\n";
    }
    else {
        die "Cannot parse infile: $infile\n";
    }
}
close $LIST;

